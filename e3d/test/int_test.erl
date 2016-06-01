%%%-------------------------------------------------------------------
%%% @author Dan Gudmundsson <dgud@erlang.org>
%%% @copyright (C) 2016, Dan Gudmundsson
%%% @doc
%%%
%%% @end
%%% Created : 20 May 2016 by Dan Gudmundsson <dgud@erlang.org>
%%%-------------------------------------------------------------------
-module(int_test).

-compile(export_all).

-define(VS, {{ 0.5,  1.0, -0.5},  %1
             { 0.5,  0.0, -0.5},  %2
             {-0.5,  0.0, -0.5},
             {-0.5,  1.0, -0.5},  %4
             {-0.5,  1.0,  0.5},
             { 0.5,  1.0,  0.5},  %6
             { 0.5,  0.0,  0.5},
             {-0.5,  0.0,  0.5}}).%8

-define(FS, [[0,1,2,3],
	     [2,7,4,3],
	     [0,5,6,1],
	     [5,4,7,6],
	     [5,0,3,4],
	     [6,7,2,1]]).


go() ->  start().

start() ->
    Ref = make_model(0, {0.0, 0.0, 0.0}),
    %% M2 = e3d_bvh:init(make_model({0.0, 1.5, 0.0})),
    test(Ref, {0.0, 1.5, 0.0}),
    test(Ref, {0.25, 0.9, 0.25}),
    ok.

test(Ref, Trans) ->
    M = make_model(1, Trans),
    Hits0 = e3d_bvh:intersect(e3d_bvh:init(Ref), e3d_bvh:init(M)),
    Hits = [I || I = #{p1:=P1,p2:=P2} <- Hits0, P1 =/= P2],
    io:format("~p: ~p~n", [length(Hits0),Hits]), %length(Hits)]),
    test2(Hits, Ref, M).

test2([], _, _) -> ok;
test2(Hits0, Ref, Other) ->
    [io:format(" ~p:~.2w  ~p:~.2w ~.5w ~s ~s~n", [M1,F1,M2,F2,P1==P2,f(P1),f(P2)])
     || #{mf1:={M1,F1},mf2:={M2,F2},p1:=P1,p2:=P2} <- Hits0],
    %Hits1 = filter_tess_edges(Hits0, Ref, Other),
    
    ok.

%filter_tess_edge({{M1,F1},{M2,F2},P1,P2}, Ref, Other) ->
    


make_model(MeshId, Trans) ->
    Fs = lists:append([ [{V1,V2,V3},{V1,V3,V4}] || [V1,V2,V3,V4] <- ?FS]),
    TFs = list_to_tuple(Fs),
    Vs = array:from_list([e3d_vec:add(Trans, Point) || Point <- tuple_to_list(?VS)]),
    format_faces(0, Fs, Vs),
    GetFace = fun({verts,Face}) -> element(Face+1, TFs);
		 (verts) -> Vs;
		 (meshId) -> MeshId;
		 (faces) -> ?FS;
		 (triangles) -> TFs
	      end,
    [{tuple_size(TFs), GetFace}].

format_faces(I, [FVs={V1,V2,V3}|Fs], Vs) ->
    io:format("~.3w ~p => [~s,~s,~s]~n",
	      [I, FVs, f(array:get(V1,Vs)),f(array:get(V2,Vs)),f(array:get(V3,Vs))]),
    format_faces(I+1, Fs,Vs);
format_faces(_,[],_) -> ok.

f({X,Y,Z}) ->
    io_lib:format("{~5.2f, ~5.2f, ~5.2f}", [X,Y,Z]);
f(X) when is_float(X) ->
    io_lib:format("~6.3f", [X]);
f(X) when is_integer(X) ->
    io_lib:format("~p", [X]).
