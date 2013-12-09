open Util
open Cpd

type edge = {
  mutable sepset: id array;
  mutable edge_cpd: cpd;
  mutable msg_waiting: int list;
}

type node = {
  id: int;
  scope: id array;
  mutable edges: node list; (* node and whether I received a msg *)
  mutable node_cpd: cpd;
  mutable saved_cpd: cpd;
}

type scheme = MaxProduct | SumProduct

let empty_edge () = {sepset=[||]; edge_cpd=empty_cpd (); msg_waiting=[]}

type tree = node * (int * int, edge) Hashtbl.t 

let lookup_edge ((_,h) : tree) (node1, node2) =
  let key = if node1.id > node2.id then node2.id, node1.id else node1.id, node2.id
  in
  Hashtbl.find h key

(* reset the cpd and mailboxes of all edges *)
let reset_edges_all ((_,h) : tree) =
  Hashtbl.iter (fun _ edge ->
    edge.edge_cpd <- empty_cpd ();
    edge.msg_waiting <- []
  ) h

(* reset the cpd and mailboxes of all edges *)
let reset_edge_mailboxes ((_,h) : tree) =
  Hashtbl.iter (fun _ edge ->
    edge.msg_waiting <- []
  ) h

let modify_edge ((_,h) : tree) (node1, node2) f =
  let key = if node1.id > node2.id 
    then node2.id, node1.id 
    else node1.id, node2.id in
  let value = try Hashtbl.find h key 
              with Not_found -> empty_edge ()
  in
  let value' = f value in
  Hashtbl.replace h key value'

let string_of_node {id; scope; edges; node_cpd} = 
  let scope_s = String.concat ", " (str_of_id_many @: Array.to_list scope) in
  let edge_ids = List.map (fun node -> string_of_int node.id) edges in
  let edges_s = String.concat "," edge_ids in
  let cpd_s = string_of_cpd node_cpd in
  Printf.sprintf "id: %d\nscope: %s\nedges: %s\ncpd: %s\n" id scope_s edges_s cpd_s

let string_of_edge {sepset; edge_cpd} =
  "sepset: "^String.concat ", " (str_of_id_many @: Array.to_list sepset)^
  "\nedge cpd: \n"^string_of_cpd edge_cpd

let tree_fold fn init ((t, _) : tree) =
  let rec loop acc node m_last_id =
    let newacc = fn acc node in
    (* loop over the edges *)
    List.fold_left (fun acc' node' ->
      match m_last_id with
      (* don't go back to where we came from *)
      | Some last_id when node'.id = last_id -> acc'
      | _ -> loop acc' node' @: Some node.id
    ) 
    newacc
    node.edges
  in loop init t None

let string_of_tree ((_, h) as tree) =
  tree_fold (fun acc n -> acc^string_of_node n) "" tree^
  Hashtbl.fold (fun (node1,node2) edge acc_string ->
    acc_string^Printf.sprintf "edge %d-%d: %s\n" node1 node2 (string_of_edge edge)
  ) h ""

let parse_clique_tree file =
  let f = read_file file in
  let ls = lines f in
  let num_nodes = int_of_string @: hd ls in
  let ls = list_drop 1 ls in
  let node_lines = list_take num_nodes ls in
  let edge_lines = list_drop num_nodes ls in
  let node_tbl = Hashtbl.create 10 in
  let count = ref 1 in
  let node_l = ref [] in
  (* add to hashtable *)
  List.iter (fun line ->
      let scope = Array.of_list @: id_of_str_many @: r_split "," line in
      let new_node = 
        {id=(!count); scope; edges=[]; node_cpd=empty_cpd (); 
         saved_cpd=empty_cpd ()} in
      Hashtbl.add node_tbl line new_node;
      node_l := new_node::(!node_l);
      count := !count + 1
    ) node_lines;
  let tree = (list_last !node_l, Hashtbl.create 10) in
  List.iter (fun edge -> 
    match r_split " -- " edge with
    | [x; y] -> 
        begin try
          let x_node = Hashtbl.find node_tbl x in
          let y_node = Hashtbl.find node_tbl y in
          x_node.edges <- y_node::x_node.edges;
          y_node.edges <- x_node::y_node.edges;
          modify_edge tree (x_node, y_node) id_fn
        with Not_found -> failwith @: "Can't find "^x^" or "^y^" in hash table" end
    | _ -> failwith "Bad split"
  ) edge_lines;
  tree

(* find the cpds relating to a node and multiply them *)
let cpd_of_node node cpd_set =
  let c = Hashtbl.fold (fun cpd _ acc_cpd ->
    if subset cpd.vars node.scope then
      (Hashtbl.remove cpd_set cpd;
      let prod = product acc_cpd cpd in
      if prod.data = [] then
        (let s1 = string_of_cpd acc_cpd in
          let s2 = string_of_cpd cpd in
        failwith @: "whoops:\n"^s1^"\n"^s2);
      prod)
    else acc_cpd
  ) cpd_set @: empty_cpd ()
  in
  (*print_endline @: string_of_cpd c;*)
  c

let tree_fill_cpds tree cpd_list =
  let h = Hashtbl.create 10 in
  List.iter (fun cpd -> Hashtbl.add h cpd ()) cpd_list;
  tree_fold (fun _ n ->
    let cpd' = cpd_of_node n h in
    n.node_cpd <- cpd';
    ()
  ) () tree;
  if Hashtbl.length h <> 0 then failwith "tree_fill_cpds: Haven't assigned all CPDs to tree"
  else tree
  
let set_sepset tree node1 =
  List.iter (fun node2 ->
    let idxs1, _, _, _ = intersect node1.scope node2.scope in
    let sepset = take_idxs idxs1 (List.length idxs1) node1.scope in
    modify_edge tree (node1, node2) (fun edge -> edge.sepset <- sepset; edge)
  ) node1.edges

let set_tree_sepsets tree = tree_fold 
   (fun _ node -> set_sepset tree node) () tree

let find_leaves tree root = 
  tree_fold (fun acc node -> 
    match node.edges with 
    | [x] when node.id <> root.id -> node::acc
    | _   -> acc
  ) [] tree

(* send a msg from node1 to node2 *)
let send_msg ?(print_send=false) ~scheme tree node1 node2 =
  let edge = lookup_edge tree (node1, node2) in
  (* get the indices of the difference between the scope and the sepset *)
  let _, _, scope_idxs, _ = intersect node1.scope edge.sepset in
  let var_set = take_idxs scope_idxs (List.length scope_idxs) node1.scope in
  let cpd_idxs = try cpd_find_idxs_arr node1.node_cpd var_set 
                 with Not_found -> let s = 
                   Printf.sprintf "Node %d scope is %s but cpd is %s. Missing vars." node1.id 
                    (sosa node1.scope)
                    (sosa node1.node_cpd.vars) in
                   failwith s
  in
  (* debug *)
  if print_send then
    (Printf.printf "Node %d scope: %s\nEdge sepset: %s\nDiff: %s\n"
      node1.id (sosa node1.scope) (sosa edge.sepset)
      (sosa var_set);
    Printf.printf "Node %d cpd:\n%s\n" node1.id (string_of_cpd node1.node_cpd));
  print_endline "marginalize";
  let msg = match scheme with
    | SumProduct -> marginalize node1.node_cpd cpd_idxs
    | MaxProduct -> marginalize_max node1.node_cpd cpd_idxs
  in
  if print_send then
    Printf.printf "Node %d marginalized cpd:\n%s\n" node1.id (string_of_cpd msg); (* debug *)
  print_endline "divide";
  let msg = div msg edge.edge_cpd in
  if print_send then
    Printf.printf "Node %d divided msg:\n%s\n" node1.id (string_of_cpd msg); (* debug *)
  edge.edge_cpd    <- msg;
  edge.msg_waiting <- node2.id::edge.msg_waiting;
  if print_send then
    Printf.printf "Node %d pre-product cpd:\n%s\n" node2.id (string_of_cpd node2.node_cpd); (* debug *)
  print_endline "product";
  node2.node_cpd <- product node2.node_cpd msg;
  print_endline "done";
  if print_send then
    Printf.printf "Node %d post-product cpd:\n%s\n" node2.id (string_of_cpd node2.node_cpd) (* debug *)

let downstream ?print_send tree =
  let root = fst tree in
  let q = Queue.create () in
  Queue.add root q;
  let rec loop node1 =
    let acc_nomsg =
      List.fold_left (fun acc_nomsg node2 ->
        let edge = lookup_edge tree (node1, node2) in
        (* send a msg *)
        if not @: List.mem node2.id edge.msg_waiting then
          (send_msg ?print_send ~scheme:SumProduct tree node1 node2;
           Queue.add node2 q);
        (* make sure we received all msgs *)
        let my_msg = List.mem node1.id edge.msg_waiting in
        if not my_msg then acc_nomsg + 1 else acc_nomsg
        ) 0 node1.edges
    in
    if acc_nomsg > 0 then 
      failwith @: "Node "^soi node1.id^"didn't receive all messages!" else
    try loop @: Queue.pop q with Queue.Empty -> ()
  in
  loop @: Queue.pop q

(* for maxproduct, we only need an upstream, and our root can vary *)
let upstream ?print_send ~scheme tree root =
  let leaves = find_leaves tree root in
  let q = Queue.create () in
  List.iter (fun n -> Queue.add n q) leaves;
  let rec loop node1 =
    let dest_node, acc_nomsg, already_sent =
      List.fold_left (fun (acc_node, acc_nomsg, acc_sent) node2 ->
          let edge = lookup_edge tree (node1, node2) in
          let msg = List.mem node1.id edge.msg_waiting in
          let already_sent = List.mem node2.id edge.msg_waiting in
          if not msg then (node2, acc_nomsg + 1, already_sent)
          else acc_node, acc_nomsg, acc_sent) 
        (node1, 0, false) 
        node1.edges
    in
    if acc_nomsg = 1 && not already_sent then 
      (send_msg ?print_send ~scheme tree node1 dest_node;
      (* don't add root so we never execute in upstream *)
      if dest_node.id <> root.id then Queue.add dest_node q else ());
    try loop @: Queue.pop q with Queue.Empty -> ()
  in
  loop @: Queue.pop q

let save_node_cpds tree =
  tree_fold (fun _ node ->
    node.saved_cpd <- node.node_cpd; ()
  ) () tree

let restore_node_cpds tree =
  tree_fold (fun _ node ->
    node.node_cpd <- node.saved_cpd; ()
  ) () tree

(* a maxproduct query root should have any one of the query vars
 * to marginalize over *)
exception Answer of node

let maxproduct_find_root tree query_vars =
  try 
    tree_fold (fun _ node ->
        let common_idxs, _, _, _ = 
          intersect node.node_cpd.vars query_vars in
        if common_idxs <> [] then
          raise @: Answer(node)
        else None)
      None
      tree
  with Answer(node) -> Some node

let normalize_tree tree =
  tree_fold (fun () node ->
    let cpd = normalize_and_real node.node_cpd in
    let cpd = cpd_to_log cpd in
    node.node_cpd <- cpd;
    ())
  () tree

