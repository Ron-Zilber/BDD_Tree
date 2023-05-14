%%% --------------------------------------------------------------------------------------------------------------------------------- %%
-module(BDD_Tree).
-author("Ron Zilber").
-export([exp_to_bdd/3, solve_bdd/2, listOfLeaves/1, reverseIteration/1]).
-define(INIT_TREE_SIZE, 0).
-define(INIT_TREE_HEIGHT, -1).
-define(INIT_NUM_OF_LEAFS, 1).
-define(MICRO_TO_MILLI, 1000).
-record(bdd_tree, {
    a_root, left_subtree, right_subtree, boolean_function, num_of_leafs, num_of_nodes, tree_height
}).
%%% --------------------------------------------------------------------------------------------------------------------------------- %%
% -------------------------- exp_to_bdd(BoolFunc, Ordering, DataStructureType) ------------------------------------------------------ %%
% Return a Bdd diagram
exp_to_bdd(BoolFunc, Ordering, DataStructureType) ->
    Start_time = erlang:timestamp(),
    % Create a tree for each possible Shannon expansion
    Trees_list = boolean_function_to_bdd_forrest(BoolFunc, DataStructureType),
    Reduced_trees_list = [reduce_bdd_tree(Tree, DataStructureType) || Tree <- Trees_list],
    Reduced_with_sizes = [get_size(Tree, DataStructureType) || Tree <- Reduced_trees_list],
    Reduced_with_sizes_and_height = [get_height(Tree, DataStructureType) || Tree <- Reduced_with_sizes],
    Reduced_with_sizes_heights_leaves_amount = [get_leaves_amount(Tree, DataStructureType)|| Tree <- Reduced_with_sizes_and_height],
    Stop_time = erlang:timestamp(),
    Elapsed_time = timer:now_diff(Stop_time, Start_time)/?MICRO_TO_MILLI,                                              %  Devide by 1000
    io:format("Function: exp_to_bdd | Data Structure: ~p | Ordering: ~p | Execution time: ~p msec~n~n", [DataStructureType, Ordering, Elapsed_time]),
    _Efficient_tree = get_efficient_tree(Reduced_with_sizes_heights_leaves_amount, DataStructureType, Ordering).

% -------------------------- solve_bdd(BddTree, Variable_Value_Pairs) --------------------------------------------------------------- %%
% Compute a boolean result for a given bdd tree and a list of values
solve_bdd(BddTree, Variable_value_pairs) ->
    Start_time = erlang:timestamp(),                                                       %  Measure the execution time of the function
    Is_map = is_map(BddTree),
    Is_record = is_record(BddTree, bdd_tree),
    if
        Is_map ->
            BoolFunc = maps:get(boolean_function, BddTree);
        Is_record ->
            #bdd_tree{boolean_function = BoolFunc} = BddTree;
        true ->
            BoolFunc = undefined
    end,
    
    Variables_list = [Variable || {Variable, _Value} <- Variable_value_pairs],                    %  Create a list of the relevant nodes
    Values_list = [Value || {_Variable, Value} <- Variable_value_pairs],                          %    Create a list of the given values
    Set_values = set_values(BoolFunc, Values_list, Variables_list),                    % Replace the variable name with the given values
    Result = logic_resolver(Set_values),                                                                    % Compute the boolean result 
    Stop_time = erlang:timestamp(),
    Elapsed_time = timer:now_diff(Stop_time, Start_time)/?MICRO_TO_MILLI,                                              %  Devide by 1000
    io:format("Function: solve_bdd | Variable - Value pairs: ~p | Result = ~p | Execution time: ~p msec~n~n", [Variable_value_pairs, Result, Elapsed_time]),
    Result.

% -------------------------- listOfLeaves(BddTree) ---------------------------------------------------------------------------------- %%
% Generate a list of path vectors from the root to each of the given tree's leaves
listOfLeaves(BddTree) ->
    Start_time = erlang:timestamp(),                                                       %  Measure the execution time of the function
    Is_map = is_map(BddTree),                               
    Is_record = is_record(BddTree, bdd_tree),
    if
        Is_map ->
            DataStructureType = map,                                                                         %  Only needed for printing
            List_of_leaves = listOfLeaves(BddTree, [{root, maps:get(a_root, BddTree)}], map),
            Return_value = tuplelize_paths(List_of_leaves, [], []),
            Stop_time = erlang:timestamp(),
            Elapsed_time = timer:now_diff(Stop_time, Start_time)/?MICRO_TO_MILLI,                                      %  Devide by 1000
            io:format("Function: listOfLeaves | Data Structure: ~p | Execution time: ~p msec~n~n", [DataStructureType, Elapsed_time]),
            Return_value;

        Is_record ->
            DataStructureType = record,
            #bdd_tree{a_root = Root} = BddTree,
            List_of_leaves = listOfLeaves(BddTree, [{root, Root}], record),
            Return_value = tuplelize_paths(List_of_leaves, [], []),
            Stop_time = erlang:timestamp(),
            Elapsed_time = timer:now_diff(Stop_time, Start_time)/?MICRO_TO_MILLI,                                      %  Devide by 1000
            io:format("Function: listOfLeaves | Data Structure: ~p | Execution time: ~p msec~n~n", [DataStructureType, Elapsed_time]),
            Return_value;
        true ->
            undefined
    end.

listOfLeaves(BddTree, Directions_list, map) ->
    if
        is_map(BddTree) ->
            Left_subtree = maps:get(left_subtree, BddTree),
            Right_subtree = maps:get(right_subtree, BddTree),

            if
                is_map(Left_subtree) ->
                    Left_list = Directions_list ++ [{left, maps:get(a_root, Left_subtree)}],
                    Left_list_next = listOfLeaves(Left_subtree, Left_list, map);
                (Left_subtree == false) or (Left_subtree == true) or (Left_subtree == 0) or
                    (Left_subtree == 1) ->
                    Left_list = Directions_list ++ [{left, Left_subtree}],
                    Left_list_next = Left_list
            end,

            if
                is_map(Right_subtree) ->
                    Right_list = Directions_list ++ [{right, maps:get(a_root, Right_subtree)}],
                    Right_list_next = listOfLeaves(Right_subtree, Right_list, map);
                (Right_subtree == false) or (Right_subtree == true) or (Right_subtree == 0) or
                    (Right_subtree == 1) ->
                    Right_list = Directions_list ++ [{right, Right_subtree}],
                    Right_list_next = Right_list
            end
    end,
    _List_of_leaves = Left_list_next ++ Right_list_next;
listOfLeaves(BddTree, Directions_list, record) ->
    if
        is_record(BddTree, bdd_tree) ->
            #bdd_tree{left_subtree = Left_subtree, right_subtree = Right_subtree} = BddTree,

            if
                is_record(Left_subtree, bdd_tree) ->
                    #bdd_tree{a_root = Left_Root} = Left_subtree,
                    Left_list = Directions_list ++ [{left, Left_Root}],
                    Left_list_next = listOfLeaves(Left_subtree, Left_list, record);
                (Left_subtree == false) or (Left_subtree == true) or (Left_subtree == 0) or
                    (Left_subtree == 1) ->
                    Left_list = Directions_list ++ [{left, Left_subtree}],
                    Left_list_next = Left_list
            end,

            if
                is_record(Right_subtree, bdd_tree) ->
                    #bdd_tree{a_root = Right_Root} = Right_subtree,
                    Right_list = Directions_list ++ [{right, Right_Root}],
                    Right_list_next = listOfLeaves(Right_subtree, Right_list, record);
                (Right_subtree == false) or (Right_subtree == true) or (Right_subtree == 0) or
                    (Right_subtree == 1) ->
                    Right_list = Directions_list ++ [{right, Right_subtree}],
                    Right_list_next = Right_list
            end
    end,
    _List_of_leaves = Left_list_next ++ Right_list_next.

% -------------------------- reverseIteration(LeafPtr) ------------------------------------------------------------------------------ %%
%  Return a list of the nodes on the shortest path from a leaf to it's tree's root
reverseIteration(LeafPtr) ->
    Start_time = erlang:timestamp(),
    Reverse_pointer = lists:reverse(LeafPtr),                                            %  Set the path nodes list to go from leaf to root
    Nodes_list = [Node || {_Atom, Node} <- Reverse_pointer],
    Stop_time = erlang:timestamp(),
    Elapsed_time = timer:now_diff(Stop_time, Start_time)/?MICRO_TO_MILLI,                                              %  Devide by 1000
    io:format("Function: reverseIteration | Execution time: ~p msec~n~n", [Elapsed_time]),
    Nodes_list.
% -------------------------- tuplelize_paths(List, Curr_list, Output_list) ---------------------------------------------------------- %%
% Get one list of various paths and seperate it to list of path tuples where each path is from the root to a leath
tuplelize_paths([], _Curr_list, Output_list) ->
    Output_list;
tuplelize_paths([H | T], Curr_list, Output_list) ->
    {Atom, Variable} = H,
    if
        Atom == root ->
            New_curr_list = [H],

            tuplelize_paths(T, New_curr_list, Output_list);
        (Variable == false) or (Variable == true) ->
            New_curr_list = Curr_list ++ [H],
            New_output_list = Output_list ++ [New_curr_list],

            tuplelize_paths(T, [], New_output_list);
        (Atom == left) or (Atom == right) ->
            New_curr_list = Curr_list ++ [H],
            tuplelize_paths(T, New_curr_list, Output_list);
        true ->
            undefined
    end.

% -------------------------- extract_vars(BoolFunc) --------------------------------------------------------------------------------- %%
% extract a list of boolean variables from a given boolean function
extract_vars(BoolFunc) ->
    ListedBoolFunc = tuplesSolver(BoolFunc),
    FlattendList = flatten(ListedBoolFunc),
    FilteredList = lists:filter(
        fun(X) -> X /= 'and' andalso X /= 'or' andalso X /= 'not' end, FlattendList
    ),
    %  Remove duplicatations of variables
    _ValuesList = lists:usort(FilteredList).

% -------------------------- tuplesSolver(Nested_tuple, Output_list) ---------------------------------------------------------------- %%
% Traverse a nested structure of tuples and lists and return a flattened list of elements
tuplesSolver([]) -> [];
tuplesSolver(Tuple) -> tuplesSolver(tuple_to_list(Tuple), []).
tuplesSolver([], List) ->
    List;
tuplesSolver([H | T], List) ->
    if
        % H is variable or operation
        not (is_tuple(H)) and not (is_list(H)) -> tuplesSolver(T, [H] ++ List);
        is_tuple(H) -> tuplesSolver(tuple_to_list(H) ++ [T], List);
        is_list(H) -> tuplesSolver(H, List);
        true -> H
    end.
% -------------------------- flatten(List) ------------------------------------------------------------------------------------------ %%

%  Returns a flat list out fo the elements in an unflvarListToTreeRecordattened input list
flatten([]) ->
    [];
flatten([H | T]) ->
    if
        %  if the head is not a list (but an element),
        (is_atom(H) or is_number(H)) -> [H] ++ flatten(T);
        %  add it to the output list and go on

        %  if the head is a list, recursively flatten it
        length(H) >= 0 -> flatten(H) ++ flatten(T);
        %  using flatten function

        %  just for completeness. Unreachable code zone (in our usage)
        true -> []
    end;
% Catches a number or an atom
flatten(Elem) ->
    [Elem].
% -------------------------- permutate_one_element(New_number, List) ---------------------------------------------------------------- %%

%  Return a list of all possible lists after adding the given element in between the elements
permutate_one_element(New_number, []) ->
    [New_number];
%  of the given list (or before / after)
permutate_one_element(New_number, List) ->
    if
        length(List) == 1 ->
            [Elem] = List,
            Output_list = [[New_number, Elem], [Elem, New_number]];
        true ->
            Indexes = lists:seq(1, length(List) + 1),
            Output_list = [expend_list(New_number, Index, List) || Index <- Indexes]
    end,
    Output_list.
% -------------------------- expend_list(Number, Index, List) ----------------------------------------------------------------------- %%

%  Returned a list with the added number in the desired index
expend_list(Number, Index, List) ->
    %  This 'if' block is used to prevent from trying to access illegal index (smaller that 1 or bigger than length(List)
    if
        (Index == 1) ->
            Returned_list = [Number] ++ List;
        % It's the N^2 -th node the (last one)
        (Index == length(List) + 1) ->
            Returned_list = List ++ [Number];
        % Just some node at the middle (index respects:     1 < index < N^2 )
        true ->
            Returned_list =
                lists:sublist(List, Index - 1) ++ [Number] ++ lists:nthtail(Index - 1, List)
    end,
    Returned_list.
% -------------------------- permutations_of_list(List) ----------------------------------------------------------------------------- %%

%  Create a lists of lists, such that each lists is a permuatation of the elements in the
permutations_of_list(List) -> permutations_of_list(List, [], 0).
%  input list
permutations_of_list([], Output_list, _Iteration) ->
    Output_list;
permutations_of_list(List, Output_list, Iteration) ->
    [H | T] = List,
    if
        Iteration < 2 ->
            New_output_list = permutate_one_element(H, Output_list),
            permutations_of_list(T, New_output_list, Iteration + 1);
        true ->
            Flattend_list = lists:flatten(Output_list),
            Splitted_list = split_list(Flattend_list, Iteration),
            Permutations_list = [
                permutate_one_element(H, Inner_list)
             || Inner_list <- Splitted_list
            ],
            Rearranged_list = split_list(lists:flatten(Permutations_list), Iteration + 1),
            permutations_of_list(T, Rearranged_list, Iteration + 1)
    end.

% -------------------------- split_list(List, K) ------------------------------------------------------------------------------------ %%
%  Split a list of length N = m * k to a list of m lists of size k
split_list(List, K) ->
    N = length(List),
    M = round(N / K),
    [lists:sublist(List, I * K + 1, K) || I <- lists:seq(0, M - 1)].

% -------------------------- boolean_function_to_bdd_forrest(BoolFunc, DataStructureType) ------------------------------------------- %%
%  Create all n! possible bdd trees from a given boolean function
boolean_function_to_bdd_forrest(BoolFunc, DataStructureType) ->
    List = extract_vars(BoolFunc),
    Permutatetions_list = lists:sort(permutations_of_list(List)),
    Trees_list = [
        build_bdd_tree(Vars_list, BoolFunc, Vars_list, [], DataStructureType)
     || Vars_list <- Permutatetions_list
    ],
    Trees_list.

% -------------------------- build_bdd_tree(List, BoolFunc, Values_list, DataStructureType) ----------------------------------------- %%
%  Build a tree (map or record) from variables and a boolean function
build_bdd_tree([H | T], BoolFunc, Variables_list, Values_list, DataStructureType) ->
    if
        %  H is the last element which means it's the last variable,
        T == [] ->
            %  his left child and right child will be the results
            Left_list = Values_list ++ [false],
            Right_list = Values_list ++ [true],
            Left_BoolFunc = set_values(BoolFunc, Left_list, Variables_list),
            Right_BoolFunc = set_values(BoolFunc, Right_list, Variables_list),
            Left_result = logic_resolver(Left_BoolFunc),
            Right_result = logic_resolver(Right_BoolFunc),

            case DataStructureType of
                map ->
                    _Results_tree = #{
                        a_root => H, left_subtree => Left_result, right_subtree => Right_result
                    };
                record ->
                    _Results_tree = #bdd_tree{
                        a_root = H, left_subtree = Left_result, right_subtree = Right_result
                    }
            end;
        %  H is not the last element, more nodes are to be added
        true ->
            Left_list = Values_list ++ [false],
            Right_list = Values_list ++ [true],

            case DataStructureType of
                map ->
                    Variables_Tree = #{a_root => H, boolean_function => BoolFunc},
                    Variables_Tree#{
                        left_subtree => build_bdd_tree(
                            T, BoolFunc, Variables_list, Left_list, DataStructureType
                        ),
                        right_subtree => build_bdd_tree(
                            T, BoolFunc, Variables_list, Right_list, DataStructureType
                        )
                    };
                record ->
                    Variables_Tree = #bdd_tree{a_root = H, boolean_function = BoolFunc},
                    Variables_Tree#bdd_tree{
                        left_subtree = build_bdd_tree(
                            T, BoolFunc, Variables_list, Left_list, DataStructureType
                        ),
                        right_subtree = build_bdd_tree(
                            T, BoolFunc, Variables_list, Right_list, DataStructureType
                        )
                    }
            end
    end.

% -------------------------- reduce_bdd_tree(Tree, DataStructureType) --------------------------------------------------------------- %%
%  Reduce a bdd_tree by the bdd reduction rules
reduce_bdd_tree(Tree, DataStructureType) ->
    case DataStructureType of
        map ->
            if
                % It is not a leaf
                is_map(Tree) ->
                    Left = reduce_bdd_tree(maps:get(left_subtree, Tree), DataStructureType),
                    Right = reduce_bdd_tree(maps:get(right_subtree, Tree), DataStructureType),

                    if
                        (Left == 0) or (Left == false) or (Left == 1) or (Left == true) ->
                            Update_left_side = Tree#{left_subtree := to_bool(Left)};
                        true ->
                            Update_left_side = Tree
                    end,

                    if
                        (Right == 0) or (Right == false) or (Right == 1) or (Right == true) ->
                            Update_right_side = Update_left_side#{right_subtree := to_bool(Right)};
                        true ->
                            Update_right_side = Update_left_side
                    end,

                    Updated_tree = Update_right_side,
                    New_left = reduce_bdd_tree(
                        maps:get(left_subtree, Updated_tree), DataStructureType
                    ),
                    New_right = reduce_bdd_tree(
                        maps:get(right_subtree, Updated_tree), DataStructureType
                    ),
                    Cond1 = is_boolean(New_left),
                    Cond2 = is_boolean(New_right),
                    Cond3 = (New_left == New_right),
                    if
                        (Cond1 and Cond2 and Cond3) ->
                            To_return_tree = New_left,
                            To_return_tree;
                        true ->
                            To_return_tree = Updated_tree,
                            To_return_tree
                    end;
                true ->
                    To_return_tree = Tree,
                    To_return_tree
            end;
        record ->
            if
                % It is not a leaf
                is_record(Tree, bdd_tree) ->
                    #bdd_tree{left_subtree = Left_subtree, right_subtree = Right_subtree} = Tree,
                    Left = reduce_bdd_tree(Left_subtree, DataStructureType),
                    Right = reduce_bdd_tree(Right_subtree, DataStructureType),
                    if
                        (Left == 0) or (Left == false) or (Left == 1) or (Left == true) ->
                            Update_left_side = Tree#bdd_tree{left_subtree = to_bool(Left)};
                        true ->
                            Update_left_side = Tree
                    end,

                    if
                        (Right == 0) or (Right == false) or (Right == 1) or (Right == true) ->
                            Update_right_side = Update_left_side#bdd_tree{
                                right_subtree = to_bool(Right)
                            };
                        true ->
                            Update_right_side = Update_left_side
                    end,

                    Updated_tree = Update_right_side,
                    #bdd_tree{
                        left_subtree = Updated_tree_left_subtree,
                        right_subtree = Updated_tree_right_subtree
                    } = Updated_tree,
                    New_left = reduce_bdd_tree(Updated_tree_left_subtree, DataStructureType),
                    New_right = reduce_bdd_tree(Updated_tree_right_subtree, DataStructureType),
                    Cond1 = is_boolean(New_left),
                    Cond2 = is_boolean(New_right),
                    Cond3 = (New_left == New_right),
                    if
                        (Cond1 and Cond2 and Cond3) ->
                            To_return_tree = New_left,
                            To_return_tree;
                        true ->
                            To_return_tree = Updated_tree,
                            To_return_tree
                    end;
                true ->
                    To_return_tree = Tree,
                    To_return_tree
            end
    end.

% -------------------------- get_size(Tree, DataStructureType) ---------------------------------------------------------------------- %%
%  Compute the size (number of nodes) of a tree and add the result to the root map
get_size(Tree, DataStructureType) ->
    case DataStructureType of
        map ->
            if
                %  It is not a leaf
                is_map(Tree) ->
                    Left = get_size(maps:get(left_subtree, Tree), DataStructureType),
                    Right = get_size(maps:get(right_subtree, Tree), DataStructureType),
                    Left_size = to_num(maps:get(num_of_nodes, Left)),
                    Right_size = to_num(maps:get(num_of_nodes, Right)),
                    _Set_size = Tree#{num_of_nodes => Left_size + Right_size + 1};
                true ->
                    _Set_size = #{num_of_nodes => ?INIT_TREE_SIZE}
            end;
        record ->
            if
                %  It is not a leaf
                is_record(Tree, bdd_tree) ->
                    #bdd_tree{left_subtree = Left_subtree, right_subtree = Right_subtree} = Tree,

                    Left = get_size(Left_subtree, DataStructureType),
                    Right = get_size(Right_subtree, DataStructureType),

                    #bdd_tree{num_of_nodes = Left_num_of_nodes} = Left,
                    #bdd_tree{num_of_nodes = Right_num_of_nodes} = Right,

                    Left_size = to_num(Left_num_of_nodes),
                    Right_size = to_num(Right_num_of_nodes),

                    _Set_size = Tree#bdd_tree{num_of_nodes = Left_size + Right_size + 1};
                true ->
                    _Set_size = #bdd_tree{num_of_nodes = ?INIT_TREE_SIZE}
            end
    end.
% -------------------------- get_size(height, DataStructureType) -------------------------------------------------------------------- %%

%  Compute the height of a tree and add the result to the root map
get_height(Tree, DataStructureType) ->
    case DataStructureType of
        map ->
            if
                % It is not a leaf
                is_map(Tree) ->
                    Left = get_height(maps:get(left_subtree, Tree), DataStructureType),
                    Right = get_height(maps:get(right_subtree, Tree), DataStructureType),

                    Left_height = to_num(maps:get(tree_height, Left)),
                    Right_height = to_num(maps:get(tree_height, Right)),

                    _Set_height = Tree#{tree_height => max(Left_height, Right_height) + 1};
                true ->
                    _Set_height = #{tree_height => ?INIT_TREE_HEIGHT}
            end;
        record ->
            if
                % It is not a leaf
                is_record(Tree, bdd_tree) ->
                    #bdd_tree{left_subtree = Left_subtree, right_subtree = Right_subtree} = Tree,

                    Left = get_height(Left_subtree, DataStructureType),
                    Right = get_height(Right_subtree, DataStructureType),

                    #bdd_tree{tree_height = Left_height} = Left,
                    #bdd_tree{tree_height = Right_height} = Right,

                    Left_height_num = to_num(Left_height),
                    Right_height_num = to_num(Right_height),

                    _Set_height = Tree#bdd_tree{
                        tree_height = max(Left_height_num, Right_height_num) + 1
                    };
                true ->
                    _Set_height = #bdd_tree{tree_height = ?INIT_TREE_HEIGHT}
            end
    end.
% -------------------------- get_leaves_amount(height, DataStructureType) ----------------------------------------------------------- %%

% Compute the number of nodes of a tree and add the result to the root map
get_leaves_amount(Tree, DataStructureType) ->
    case DataStructureType of
        map ->
            if
                % It is not a leaf
                is_map(Tree) ->
                    Left = get_leaves_amount(maps:get(left_subtree, Tree), DataStructureType),
                    Right = get_leaves_amount(maps:get(right_subtree, Tree), DataStructureType),
                    Left_num_of_leaves = to_num(maps:get(num_of_leafs, Left)),
                    Right_num_of_leaves = to_num(maps:get(num_of_leafs, Right)),
                    _Set_num_of_leaves = Tree#{
                        num_of_leafs => (Left_num_of_leaves + Right_num_of_leaves)
                    };
                true ->
                    _Set_num_of_leaves = #{num_of_leafs => ?INIT_NUM_OF_LEAFS}
            end;
        record ->
            if
                % It is not a leaf
                is_record(Tree, bdd_tree) ->
                    #bdd_tree{left_subtree = Left_subtree, right_subtree = Right_subtree} = Tree,
                    Left = get_leaves_amount(Left_subtree, DataStructureType),
                    Right = get_leaves_amount(Right_subtree, DataStructureType),

                    #bdd_tree{num_of_leafs = Left_num_of_leaves} = Left,
                    #bdd_tree{num_of_leafs = Right_num_of_leaves} = Right,

                    Left_num_of_leaves_num = to_num(Left_num_of_leaves),
                    Right_num_of_leaves_num = to_num(Right_num_of_leaves),

                    _Set_num_of_leaves = Tree#bdd_tree{
                        num_of_leafs = (Left_num_of_leaves_num + Right_num_of_leaves_num)
                    };
                true ->
                    _Set_num_of_leaves = #bdd_tree{num_of_leafs = ?INIT_NUM_OF_LEAFS}
            end
    end.
% -------------------------- get_efficient_tree(Trees_list , DataStructureType, Ordering) ------------------------------------------- %%

%  Return the most efficient tree by a given criteria
get_efficient_tree(Trees_list, DataStructureType, Ordering) ->
    case Ordering of
        num_of_nodes ->
            case DataStructureType of
                map ->
                    Values_list = [maps:get(num_of_nodes, Tree) || Tree <- Trees_list];
                record ->
                    Values_list = [
                        Num_of_nodes
                     || #bdd_tree{num_of_nodes = Num_of_nodes} <- Trees_list
                    ]
            end;
        num_of_leafs ->
            case DataStructureType of
                map ->
                    Values_list = [maps:get(num_of_leafs, Tree) || Tree <- Trees_list];
                record ->
                    Values_list = [
                        Num_of_leaves
                     || #bdd_tree{num_of_leafs = Num_of_leaves} <- Trees_list
                    ]
            end;
        tree_height ->
            case DataStructureType of
                map ->
                    Values_list = [maps:get(tree_height, Tree) || Tree <- Trees_list];
                record ->
                    Values_list = [
                        Tree_height
                     || #bdd_tree{tree_height = Tree_height} <- Trees_list
                    ]
            end;
        _ ->
            Values_list = undefined,
            undefined
    end,

    Minimal_value = lists:min(Values_list),
    Index_of_most_efficient_tree = get_index(Values_list, Minimal_value),
    _Returned_Tree = lists:nth(Index_of_most_efficient_tree, Trees_list).

% -------------------------- get_index(List, Num) ----------------------------------------------------------------------------------- %%

% Return the index of a given element in a given list
get_index(List, Num) -> get_index(List, Num, 1).
get_index([], _Num, _Index) -> undefined;
get_index([Num | _T], Num, Index) -> Index;
get_index([_H | T], Num, Index) -> get_index(T, Num, Index + 1).

% -------------------------- to_bool(value) ----------------------------------------------------------------------------------------- %%
%  return a boolean value in the form of false/true
to_bool(0) -> false;
to_bool(false) -> false;
to_bool(1) -> true;
to_bool(true) -> true.
% -------------------------- to_num(value) ------------------------------------------------------------------------------------------ %%
to_num(false) -> ?INIT_TREE_HEIGHT;
to_num(true) -> ?INIT_TREE_HEIGHT;
to_num(Num) -> Num.
% -------------------------- logic_resolver(BooleanFunc) ---------------------------------------------------------------------------- %%
% Compute the boolean value of a boolean expression with 0/1, false/true values
logic_resolver({'and', {X1, X2}}) ->
    logic_resolver(X1) and logic_resolver(X2);
logic_resolver({'or', {X1, X2}}) ->
    logic_resolver(X1) or logic_resolver(X2);
logic_resolver({'not', X}) ->
    not (logic_resolver(X));
logic_resolver(X) ->
    if
        ((X == 1) or (X == true)) ->
            true;
        ((X == 0) or (X == false)) ->
            false;
        true ->
            undefined
    end.

%%% ------------------------ set_values(BoolFunc, Values_list) ---------------------------------------------------------------------- %%

% replace atoms x1, x2, ,,, xn , etc to '0' or '1' values form the list
set_values(BoolFunc, Values_list, Variables_list) ->
    Indexes = lists:seq(1, length(Values_list)),
    set_values(BoolFunc, Variables_list, Values_list, Indexes).
set_values(BoolFunc, _Variables_list, _Values_list, []) ->
    BoolFunc;
set_values(BoolFunc, Variables_list, Values_list, [H | T]) ->
    New_BoolFunc = set_value(BoolFunc, Variables_list, Values_list, H),
    set_values(New_BoolFunc, Variables_list, Values_list, T);
set_values(BoolFunc, Variables_list, Values_list, H) ->
    New_BoolFunc = set_value(BoolFunc, Variables_list, Values_list, H),
    set_values(New_BoolFunc, Variables_list, Values_list, []).

%%% ------------------------ set_value(BoolFunc, Variables_list, Values_list, Index) ------------------------------------------------ %%
% Replace a variable name in a boolean expression to the corresponding boolean value
set_value(BoolFunc, Variables_list, Values_list, Index) ->
    String_bool_func = lists:flatten(io_lib:format("~p", [BoolFunc])),
    With_dot = string:join([String_bool_func, "."], ""),
    Variable_Name = io_lib:format("~p", [(lists:nth(Index, Variables_list))]),
    [New_Value] = io_lib:format("~p", [(lists:nth(Index, Values_list))]),
    Unicode_data = string:replace(With_dot, Variable_Name, New_Value, all),
    Str = unicode:characters_to_list(Unicode_data, {utf16, little}),
    {ok, Ts, _} = erl_scan:string(Str),
    {ok, Tuple} = erl_parse:parse_term(Ts),
    Tuple.
%%% --------------------------------------------------------------------------------------------------------------------------------- %%

%%% --------------------------------------------------------------------------------------------------------------------------------- %%
%%% ------------------------ Thank you for reading ---------------------------------------------------------------------------------- %%
