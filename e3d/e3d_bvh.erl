%%%-------------------------------------------------------------------
%%% @author Dan Gudmundsson <dgud@erlang.org>
%%% @copyright (C) 2016, Dan Gudmundsson
%%% @doc Bounding Volume Hierarchy
%%%      Simple tree for ray-testing or intersections
%%%      Output can be term() tree or a binary() for OpenCL or for parallel raytracing
%%% @end
%%% Created : 29 Apr 2016 by Dan Gudmundsson <dgud@erlang.org>
%%%-------------------------------------------------------------------
-module(e3d_bvh).
-export([init/1, init/2,
	 ray/2, ray/4, ray_trace/2,
	 intersect/2,
	 hit_triangle/5, tri_intersect/4]).

-include("e3d.hrl").

-type hit() ::
	#{t  => float(),
	  b1 => float(),
	  b2 => float(),
	  mesh => integer(),
	  face => integer()}.

-type tri_intersect() ::  %% The (new) edge that intersects the two triangle
	#{p1    => e3d_point(),  %% Edge point 1
	  p2    => e3d_point(),  %% Edge point 2
	  mf1   => {Mesh::integer(),Face::integer()},  %% Mesh and Face of P1
	  mf2   => {Mesh::integer(),Face::integer()},  %% Mesh and Face of P2
	  other => {Mesh::integer(),Face::integer()}}. %% Intersecting Mesh and Face

-type e3d_bvh() :: #{vs => #{integer()=>array:array()} | binary(),
		     ns => tree() | binary()}.

-type e3d_compiled() :: {bvh, compiled, atom()}.

-type leaf() :: #{bb    => e3d_bbox(),
		  c2    => e3d_point(),  %% CenterPos*2
		  vs    => {integer(), integer(), integer()}, %% Vertex indices
		  mesh  => integer(),
		  index => integer()}.

-type tree_node() :: #{bb    => e3d_bbox(),
		       left  => tree_node() | leaf(),
		       right => tree_node() | leaf()}.

-type tree() :: leaf() | tree_node().

-define(F32, 32/float-native).
-define(I32, 32/signed-native).
-define(U32, 32/unsigned-native).

-define(EPSILON, 0.000001).
-define(IsLeaf(NodeData), ((NodeData) band 16#80000000) > 1).
-define(GetSkipIndex(NodeData), ((NodeData) band 16#7fffffff)).
-define(ENTRY_SZ, 8*4).

-spec init([{NoFs :: integer(), GetVs :: function()}] | e3d_compiled()) -> e3d_bvh().
init({bvh, compiled, Mod}) ->
    #{vs=>Mod:verts(), ns=>Mod:tree()};
init(Data) ->
    init(Data, []).

-spec init([{NoFs :: integer(), GetVs :: function()}], []) -> e3d_bvh() | e3d_compiled().
init(FaceData, Opts)  ->
    TT  = proplists:get_value(treetype, Opts, 4),
    Eps = proplists:get_value(epsilon, Opts, ?EPSILON),
    Nodes= lists:foldl(fun({N, GetVs}, Acc) ->
			       make_nodes(N-1, GetVs, Eps, Acc)
		       end, [], FaceData),
    Root = build_hierarchy(Nodes, TT),
    case proplists:get_value(binary, Opts, false) of
	false  ->
	    Vs0 = [{GetVs(meshId),GetVs(verts)}
		   || {_, GetVs} <- FaceData],
	    Vs = maps:from_list(Vs0),
	    case proplists:get_value(compile, Opts, false) of
		false ->
		    #{vs=>Vs, ns=>Root};
		true ->
		    compile(Root, Vs)
	    end;
	true ->
	    DeepVs0 = [GetVs(verts) || {_, GetVs} <- FaceData],
	    {VsBin,VsOffsets} = build_bin(DeepVs0, <<>>, [0]),
	    {_,A} = build_array(Root, 0, array:new(), VsOffsets),
	    #{vs=>VsBin, ns=>list_to_binary(array:to_list(A))}
    end.

%%--------------------------------------------------------------------
%% @doc Creates a ray
%% @end
%%--------------------------------------------------------------------
-spec ray(e3d_point(), e3d_vector()) -> e3d_ray().
ray(Orig, Vector) ->
    ray(Orig, Vector, ?EPSILON*10, ?E3D_INFINITY).

-spec ray(e3d_point(), e3d_vector(), float(), float()) -> e3d_ray().
ray(Orig, Vector, MinT, MaxT) ->
    #ray{o=Orig, d=Vector, n=MinT, f=MaxT}.

%%--------------------------------------------------------------------
%% @doc Cast a ray on BVH
%% @end
%%--------------------------------------------------------------------
-spec ray_trace(e3d_ray(), e3d_bvh()) -> false | hit().
ray_trace(Ray = #ray{d=Dir}, #{ns:=#{}=Root, vs:=Vs}) ->
    Swap  = e3d_bv:inv_sign(Dir),
    {_, Hit} = erl_ray_trace(Ray, false, Swap, Root, Vs),
    Hit;
ray_trace(Ray = #ray{d=Dir}, #{vs:=VsBin, ns:=Bin}) when is_binary(VsBin) ->
    Swap  = e3d_bv:inv_sign(Dir),
    {_, Hit} = bin_ray_trace(0, Ray, false, Swap, Bin, VsBin, 1),
    Hit.

%%--------------------------------------------------------------------
%% @doc Intersect two BVH's and return new edges on intersecting
%% triangles if any.
%% @end
%% --------------------------------------------------------------------
-spec intersect(e3d_bvh(), e3d_bvh()) -> [tri_intersect()].
intersect(#{ns:=N1, vs:=Vs1}, #{ns:=N2, vs:=Vs2}) ->
    intersect_1(N1, N2, Vs1, Vs2, []).

%% --------------------------------------------------------
%% Internals
%% --------------------------------------------------------

make_nodes(N, GetVs, Eps, Ns) when N >= 0 ->
    Tri  = GetVs({verts, N}),
    Mesh = GetVs(meshId),
    {V1,V2,V3} = get_tri(Tri, GetVs(verts)),
    {Min,Max} = BB = e3d_bv:box([V1,V2,V3], Eps),
    C2 = e3d_vec:add(Min,Max),
    Node = #{bb=>BB, c2=>C2, vs=>Tri, mesh=>Mesh, index=>N},
    make_nodes(N-1, GetVs, Eps, [Node|Ns]);
make_nodes(_, _, _, Acc) ->
    Acc.

build_hierarchy([Node], _TT) -> Node;
build_hierarchy(Nodes, TT) ->
    Split = find_best_split(Nodes),
    case partition(Split, Nodes) of
	{L1,L2} ->
	    #{bb:=LBB} = Left  = build_hierarchy(L1, TT),
	    #{bb:=RBB} = Right = build_hierarchy(L2, TT),
	    #{bb=>e3d_bv:union(LBB,RBB), left=>Left, right=>Right}
    end.

find_best_split(Nodes) ->
    C0 = lists:foldl(fun(#{c2:=C}, Acc) -> e3d_vec:add(C,Acc) end, e3d_vec:zero(), Nodes),
    %io:format("~p~n ~p~n",[C0, Nodes]),
    {Mx,My,Mz} = Mean2 = e3d_vec:divide(C0, float(length(Nodes))),
    {Vx,Vy,Vz} = lists:foldl(fun(#{c2:=C}, Acc) ->
				     {X,Y,Z} = e3d_vec:sub(C, Mean2),
				     V = {X*X,Y*Y,Z*Z},
				     e3d_vec:add(Acc, V)
			     end, e3d_vec:zero(), Nodes),
    %io:format("SPLIT: ~p < ~p~n", [{Vx,Vy,Vz}, Mean2]),
    if Vx =:= Vy, Vx =:= Vz -> %% Perfectly centered use mean
	    %io:format("Use Mean~n",[]),
	    if abs(Mx) < abs(Mz), abs(My) < abs(Mz) -> fun(#{c2:={_,_,V}}) -> V < Mz end;
	       abs(Mx) < abs(My) -> fun(#{c2:={_,V,_}}) -> V < My end;
	       true -> fun(#{c2:={V,_,_}}) -> V < Mx end
	    end;
       Vx < Vz, Vy =< Vz ->  %io:format("split Z ~p~n",[Mz]),
	    fun(#{c2:={_,_,V}}) -> V < Mz end;
       Vx < Vy ->           %io:format("split Y ~p~n",[My]),
	    fun(#{c2:={_,V,_}}) -> V < My end;
       true ->              %io:format("split X ~p~n",[Mx]),
	    fun(#{c2:={V,_,_}}) -> V < Mx end
    end.

partition(_, [N1,N2]) ->
    {[N1],[N2]};
partition(Split, Nodes) ->
    lists:partition(Split, Nodes).

build_array(#{bb:=BB, left:=L, right:=R}, Offset0, Array0, VsOffsets) ->
    {Offset1, Array1} = build_array(L, Offset0+1,  Array0, VsOffsets),
    {Offset2, Array2} = build_array(R, Offset1,  Array1, VsOffsets),
    {{V0x,V0y,V0z},{V1x,V1y,V1z}} = BB,
    Bin = <<V0x:?F32, V0y:?F32, V0z:?F32,V1x:?F32, V1y:?F32, V1z:?F32, Offset2:?U32, 0:?U32>>,
    {Offset2, array:set(Offset0, Bin, Array2)};
build_array(#{index:=I, mesh:=MeshId, vs:={V1,V2,V3}}, Offset0, Array0, VsOffsets) ->
    Leaf = (Offset0+1) bor 16#80000000,
    Pad = 0,
    VsOff = lists:nth(MeshId+1, VsOffsets),
    Bin = <<(V1+VsOff):?U32, (V2+VsOff):?U32, (V3+VsOff):?U32,
	    MeshId:?U32, I:?U32, Pad:?U32,  Leaf:?U32, Pad:?U32>>,
    {Offset0+1, array:set(Offset0, Bin, Array0)}.

build_bin([Top|DeepVs0], Acc0, [Prev|_]=Indx) ->
    Acc = array:foldl(fun(_, {X,Y,Z}, Acc) -> <<Acc/binary, X:?F32, Y:?F32, Z:?F32>> end,
		      Acc0, Top),
    Sz = (byte_size(Acc) - byte_size(Acc0)) div 12,
    build_bin(DeepVs0, Acc, [Prev+Sz|Indx]);
build_bin([], Acc, Indx) ->
    {Acc, lists:reverse(Indx)}.

erl_ray_trace(Ray0, Hit0, Swap, #{bb:=BB, left:=L, right:=R}, Vs) ->
    case e3d_bv:hit(Ray0, Swap, BB) of
	true ->
	    {Ray, Hit} = erl_ray_trace(Ray0, Hit0, Swap, L, Vs),
	    erl_ray_trace(Ray, Hit, Swap, R, Vs);
	false ->
	    {Ray0, Hit0}
    end;
erl_ray_trace(R, H, _Swap, #{mesh:=Mesh, index:=I, vs:=Face}, Vs) ->
    hit_triangle(R, H, get_tri(Mesh, Face, Vs), Mesh, I).

bin_ray_trace(Current, Ray0, Hit0, Swap, Bin, Vs, Level) ->
    Offset = (Current*?ENTRY_SZ),
    case Bin of
	<<_:Offset/binary, Data:(6*4)/binary, NodeData:?U32, _/binary>> ->
	    case ?IsLeaf(NodeData) of
		true  ->
		    {Ray, Hit} = bin_tri_hit(Ray0, Hit0, Data, Vs),
		    bin_ray_trace(Current+1, Ray, Hit, Swap, Bin, Vs, Level);
		false ->
		    case bin_bb_intersect(Ray0, Swap, Data) of
			true  -> bin_ray_trace(Current+1, Ray0, Hit0, Swap, Bin, Vs, Level+1);
			false -> bin_ray_trace(NodeData, Ray0, Hit0, Swap, Bin, Vs, Level+1)
		    end
	    end;
	<<_:Offset/binary>> -> {Ray0, Hit0}
    end.

bin_tri_hit(Ray, Hit, <<V1:?U32, V2:?U32, V3:?U32, MeshId:?U32, I:?U32, _:?U32>>, Vs) ->
    <<_:V1/binary-unit:96, X1:?F32, Y1:?F32, Z1:?F32, _/binary>> = Vs,
    <<_:V2/binary-unit:96, X2:?F32, Y2:?F32, Z2:?F32, _/binary>> = Vs,
    <<_:V3/binary-unit:96, X3:?F32, Y3:?F32, Z3:?F32, _/binary>> = Vs,
    V1p = {X1,Y1,Z1},
    V2p = {X2,Y2,Z2},
    V3p = {X3,Y3,Z3},
    hit_triangle(Ray, Hit, {V1p, V2p, V3p}, MeshId, I).

bin_bb_intersect(Ray, Swap, <<MIx:?F32, MIy:?F32, MIz:?F32, MAx:?F32, MAy:?F32, MAz:?F32>>) ->
    e3d_bv:hit(Ray, Swap, {{MIx,MIy,MIz},{MAx,MAy,MAz}}).

%% Woop JCGT 2(1)
%% http://jcgt.org/published/0002/01/05/paper.pdf
hit_triangle(#ray{o=Orig, d=Dir}=Ray, Hit0, % {_,_,Kz} = Order, {_,_,Sz} = Shear,
	     {TA,TB,TC}, Mesh, Face) ->
    {_,_,Kz} = Order = order(Dir),
    {_,_,Sz} = Shear = shear(Dir,Order),

    %% Calc vs relative ray origin
    A = e3d_vec:sub(TA, Orig),
    B = e3d_vec:sub(TB, Orig),
    C = e3d_vec:sub(TC, Orig),
    %% Shear and scale vs
    {Ax, Ay} = shear_scale(A, Shear, Order),
    {Bx, By} = shear_scale(B, Shear, Order),
    {Cx, Cy} = shear_scale(C, Shear, Order),
    %% Calc scaled barycentric coords
    case bary_centric(Ax,Ay,Bx,By,Cx,Cy) of
	false -> {Ray,Hit0};
	{U,V,W,Det} ->
	    %% Calc scaled z-coords and use them for hit dist
	    Az = Sz*element(Kz, A),
	    Bz = Sz*element(Kz, B),
	    Cz = Sz*element(Kz, C),
	    T = U*Az + V*Bz + W*Cz,
	    case Hit0 of
		%% ifdef backface culling
		_ when T < 0.0 ->
		    {Ray, Hit0};
		#{t:=HitT} when T > (HitT*Det) ->
		    {Ray, Hit0};
		%% non backface culling
		%% (Det > 0.0 andalso T < 0.0) orelse
		%% (Det < 0.0 andalso T > 0.0) orelse
		%% -1.0*T > -1.0*Hit0#hit.t*Det) -> Hit0;
		_ ->
		    RcpDet = 1.0 / Det,
		    Hit = #{t=>Far=T*RcpDet,
			    b1=>V*RcpDet,
			    b2=>W*RcpDet,
			    mesh=>Mesh,
			    face=>Face},
		    {Ray#ray{f=Far}, Hit}
	    end
    end.

bary_centric(Ax,Ay,Bx,By,Cx,Cy) ->
    U = Cx*By-Cy*Bx,
    V = Ax*Cy-Ay*Cx,
    W = Bx*Ay-By*Ax,
    %% Backface culling
    if U < 0.0; V < 0.0; W < 0.0 -> false;
       %% %% non backface culling
       %% if (U < 0.0 orelse V < 0.0 orelse W < 0.0) andalso
       %%    (U > 0.0 orelse V > 0.0 orelse W > 0.0) ->
       %% 	    false;
       true ->
	    case U+V+W of
		0.0 -> false;
		Det -> {U,V,W,Det}
	    end
    end.

order({X0,Y0,Z0}) ->
    %% Calc dimension where ray dir is maximal
    %% and swap the other to preserve winding
    X = abs(X0), Y = abs(Y0), Z = abs(Z0),
    if X >= Y ->
	    case X >= Z of
		true  when X0 >= 0 -> {2, 3, 1};
		true               -> {3, 2, 1};
		false when Z0 >= 0 -> {1, 2, 3};
		false              -> {2, 1, 3}
	    end;
       true ->
	    case Y >= Z of
		true  when Y0 >= 0 -> {3, 1, 2};
		true               -> {1, 3, 2};
		false when Z0 >= 0 -> {1, 2, 3};
		false              -> {2, 1, 3}
	    end
    end.

shear(Dir, {OX,OY,OZ}) ->
    DZ = element(OZ,Dir),
    {element(OX,Dir)/DZ, element(OY,Dir)/DZ, 1.0/DZ}.

shear_scale(Vec, {Sx,Sy,_Sz}, {OX,OY,OZ}) ->
    Az = element(OZ, Vec),
    {element(OX,Vec)-Sx*Az, element(OY,Vec)-Sy*Az}.

%%-------------------------------------------------------

intersect_1(#{bb:=BB1}=N1, #{bb:=BB2}=N2, Vs1, Vs2, Acc) ->
    case e3d_bv:intersect(BB1, BB2) of
	true  -> intersect_2(N1, N2, Vs1, Vs2, Acc);
        false -> Acc
    end.

intersect_2(#{left:=L1,right:=R1},  #{left:=L2,right:=R2}, Vs1, Vs2, Acc0) ->
    Acc1 = intersect_1(L1, L2, Vs1, Vs2, Acc0),
    Acc2 = intersect_1(L1, R2, Vs1, Vs2, Acc1),
    Acc3 = intersect_1(L2, R1, Vs2, Vs1, Acc2),
    intersect_1(R2, R1, Vs2, Vs1, Acc3);
intersect_2(#{left:=L1,right:=R1}, Leaf, Vs1, Vs2, Acc) ->
    intersect_1(L1, Leaf, Vs1, Vs2, intersect_1(R1, Leaf, Vs1, Vs2, Acc));
intersect_2(Leaf, #{left:=L1,right:=R1}, Vs1, Vs2, Acc) ->
    intersect_1(L1, Leaf, Vs2, Vs1, intersect_1(R1, Leaf, Vs2, Vs1, Acc));
intersect_2(#{mesh:=Mesh1, index:=I1, vs:=F1},
	    #{mesh:=Mesh2, index:=I2, vs:=F2},
	    Vs1, Vs2, Acc) ->
    T1 = get_tri(Mesh1, F1, Vs1),
    T2 = get_tri(Mesh2, F2, Vs2),
    %% io:format("~p ~p:~p:~p~n~p ~p:~p:~p",
    %% 	      [Mesh1, I1  div 2, I1, T1, Mesh2, I2 div 2, I2, T2]),
    case tri_intersect(T1, T2, {Mesh1,I1}, {Mesh2,I2}) of
	false ->
	    %% io:format(" => Miss~n~n"),
	    Acc;
	Intersect  ->
	    %% io:format(" => Hit~n~n"),
	    [Intersect|Acc]
    end.

%% Möller (realtimerendering, page 590 in 2nd edition)
tri_intersect({V0,V1,V2}, {U0,U1,U2}, F1, F2) ->
    E1 = e3d_vec:sub(V1, V0),
    E2 = e3d_vec:sub(V2, V0),
    N1 = e3d_vec:cross(E1,E2),
    D1 = -e3d_vec:dot(N1,V0),
    %% Plane equation 1: N1x+D1=0

    %% Put U0,U1,U2 in Plane eq 1 to compute signed distances to the plane
    Du0 = eps(e3d_vec:dot(N1,U0)+D1),
    Du1 = eps(e3d_vec:dot(N1,U1)+D1),
    Du2 = eps(e3d_vec:dot(N1,U2)+D1),

    Du0Du1 = Du0*Du1,
    Du0Du2 = Du0*Du2,

    case Du0Du1 > 0.0 andalso Du0Du2 > 0.0 of
	true -> false; %% No intersection occur, triangle above or under plane 1
	false ->
	    E3 = e3d_vec:sub(U1, U0),
	    E4 = e3d_vec:sub(U2, U0),
	    N2 = e3d_vec:cross(E3,E4),
	    D2 = -e3d_vec:dot(N2,U0),
	    %% Plane equation 2: N2x+D2=0
	    Dv0 = eps(e3d_vec:dot(N2,V0)+D2),
	    Dv1 = eps(e3d_vec:dot(N2,V1)+D2),
	    Dv2 = eps(e3d_vec:dot(N2,V2)+D2),

	    Dv0Dv1 = Dv0*Dv1,
	    Dv0Dv2 = Dv0*Dv2,

	    case Dv0Dv1 > 0.0 andalso Dv0Dv2 > 0.0 of
		true -> false; %% No intersection occur, triangle above or under plane 2
		false ->
		    %% Compute direction of intersection line
		    D = e3d_vec:cross(N1,N2),
		    Index = largest_dir(D),
		    %% Compute interval for triangle 1
		    case tri_intvals(V0,V1,V2, Index, Dv0, Dv1, Dv2, Dv0Dv1, Dv0Dv2) of
			true ->
			    %% case coplanar_tri_tri(N1,V0,V1,V2,U0,U1,U2) of
			    %% 	true -> {A1,A2};
			    %% 	false -> false
			    %% end;
			    %% Coplanar we don't care about those faces
			    coplanar;
			{ISect1,A1,A2} ->
			    {ISect2,B1,B2} = tri_intvals(U0,U1,U2, Index, Du0, Du1, Du2, Du0Du1, Du0Du2),
			    %% io:format("~p  ~p~n",[sort2(ISect1),sort2(ISect2)]),
			    %% io:format("A1:~s~nA2:~s~nB1:~s~nB2:~s~n", [f(A1),f(A2),f(B1),f(B2)]),
			    pick_points(sort2(ISect1), sort2(ISect2), A1, A2, B1, B2, F1, F2)
		    end
	    end
    end.

tri_intvals(V0, V1, V2, Index, D0, D1, D2, DOD1, _DOD2)
  when DOD1 > 0.0 ->
    %% here we know that D0D2<=0.0
    %%  that is D0, D1 are on the same side, D2 on the other or on the plane
    isect2(V2,V0,V1,Index,D2,D0,D1);
tri_intvals(V0, V1, V2, Index, D0, D1, D2, _DOD1, DOD2)
  when DOD2 > 0.0 ->
    %% here we know that d0d1<=0.0
    isect2(V1,V0,V2,Index,D1,D0,D2);
tri_intvals(V0, V1, V2, Index, D0, D1, D2, _DOD1, _DOD2)
  when (D1*D2>0.0) orelse  D0/=0.0 ->
    isect2(V0,V1,V2,Index,D0,D1,D2);
tri_intvals(V0, V1, V2, Index, D0, D1, D2, _DOD1, _DOD2)
  when D1/=0.0 ->
    isect2(V1,V0,V2,Index,D1,D0,D2);
tri_intvals(V0, V1, V2, Index, D0, D1, D2, _DOD1, _DOD2)
  when D2/=0.0 ->
    isect2(V2,V0,V1,Index,D2,D0,D1);
tri_intvals(_V1, _V2, _V3, _Index, _D0, _D1, _D2, _DOD1, _DOD2) ->
    %% triangles are coplanar
    true.

isect2(V0, V1, V2, Index, D0, D1, D2) ->
    Tmp0 = D0/(D0-D1),
    VV0 = element(Index, V0),
    VV1 = element(Index, V1),
    VV2 = element(Index, V2),
    Isect0 = VV0+(VV1-VV0)*Tmp0,
    Diff00 = e3d_vec:sub(V1,V0),
    Diff01 = e3d_vec:mul(Diff00, Tmp0),
    P0 = e3d_vec:add(V0, Diff01),
    Tmp1 = D0/(D0-D2),
    Isect1 = VV0+(VV2-VV0)*Tmp1,
    Diff10 = e3d_vec:sub(V2,V0),
    Diff11 = e3d_vec:mul(Diff10, Tmp1),
    P1 = e3d_vec:add(V0, Diff11),
    {{Isect0, Isect1}, P0, P1}.

pick_points({IS10,IS11,_}, {IS20,IS21,_}, _A1, _A2, _B1, _B2, _F1, _F2)
  when (IS11 < IS20) orelse (IS21 < IS10) ->
    false;
pick_points({IS10,IS11,Min1}, {IS20,IS21,Min2}, A1, A2, B1, B2, F1, F2)
  when IS20 < IS10 ->
    {MF1,P1} = {F1, pick_point(Min1,A1,A2)},
    {MF2,P2} = if IS21 < IS11 -> {F2, pick_point(Min2,B2,B1)};
		  true -> {F1, pick_point(Min1,A2,A1)}
	       end,
    #{mf1=>MF1, mf2=>MF2, p1=>P1, p2=>P2, other=>F2};
pick_points({_IS10,IS11,Min1}, {_IS20,IS21,Min2}, A1, A2, B1, B2, F1, F2) ->
    {MF1,P1} = {F2,pick_point(Min2, B1, B2)},
    {MF2,P2} = if IS21 > IS11 -> {F1,pick_point(Min1, A2, A1)};
		  true -> {F2,pick_point(Min2, B2, B1)}
	       end,
    #{mf1=>MF1, mf2=>MF2, p1=>P1, p2=>P2, other=>F1}.

largest_dir({X,Y,Z}) ->
    AX=abs(X), AY=abs(Y), AZ=abs(Z),
    if AY < AZ, AX < AZ -> 3;
       AX < AY -> 2;
       true -> 1
    end.

pick_point(true, P1, _P2) -> P1;
pick_point(false,_P1, P2) -> P2.

sort2({A,B}) when A > B ->
    {B,A,false};
sort2({A,B}) ->
    {A,B,true}.

eps(V) when abs(V) < ?EPSILON -> 0.0;
eps(V) -> V.

get_tri(MeshId, {V1,V2,V3}, AllVs) ->
    Vs = maps:get(MeshId,AllVs),
    %%{element(V1+1, Vs),element(V2+1, Vs),element(V3+1, Vs)}.
    {array:get(V1, Vs),array:get(V2, Vs),array:get(V3, Vs)}.

get_tri({V1,V2,V3}, Vs) ->
    {array:get(V1, Vs),array:get(V2, Vs),array:get(V3, Vs)}.

%%%% Compile to constants

compile(Tree, Verts) ->
    Mod = "tmp_" ++ ?MODULE_STRING ++ integer_to_list(erlang:unique_integer([positive,monotonic])),
    H1 = io_lib:format("-module(~s).~n-compile(export_all).~n~n", [Mod]),
    F1 = io_lib:format("tree() ->~n ~w.~n~n", [Tree]),
    F2 = io_lib:format("verts() ->~n ~w.~n", [Verts]),
    {Dir,File} = try
		     BD = filename:basedir(user_cache, "e3d"),
		     F = filename:join(BD, Mod),
		     ok = filelib:ensure_dir(F),
		     {BD, F}
		 catch _:_ ->
			 {"/tmp", filename:join("/tmp", Mod)}
		 end,
    ok = file:write_file(File  ++ ".erl", iolist_to_binary([H1,F1,F2])),
    {ok, Module} = c:c(File, [{outdir, Dir}]),
    ok = file:delete(File ++ ".erl"),
    ok = file:delete(File ++ ".beam"),
    {bvh, compiled, Module}.
