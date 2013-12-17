open Util
open CliqueTree
open Cpd
open Queries
open Unix

type actions = Print_Tree
             | Print_CPDs
             | Inference
             | MaxProductInference
             | ComputeJoint

type params_t = {
  mutable action: actions;
  mutable network_file: string;
  mutable cpd_file: string;
  mutable cliquetree_file: string;
  mutable queries_file: string;
  mutable debug_send: bool;
  mutable print_tree: bool;
  mutable incremental: bool;
  mutable time: bool;
}


let params = {
  action=Inference;
  network_file="";
  cpd_file="";
  cliquetree_file="";
  queries_file="";
  debug_send=false;
  print_tree=false;
  incremental=false;
  time=false;
  }

let parse_cmd_line () =
  let files = ref [] in
  let param_specs = Arg.align
   [
    "--incremental", Arg.Unit (fun () -> params.incremental <- true),
     "Incremental computation";
    "--max", Arg.Unit (fun () -> params.action <- MaxProductInference),
     "Use max-product instead of sum-product";
    "--joint", Arg.Unit (fun () -> params.action <- ComputeJoint),
     "Compute the joint for a certain assignment";
    "--print_cpds", Arg.Unit (fun () -> params.action <- Print_CPDs),
     "Only print the CPDs";
    "--print_init_tree", Arg.Unit (fun () -> params.action <- Print_Tree),
     "Only print the initialized clique tree";
    "--print_tree", Arg.Unit (fun () -> params.print_tree <- true),
     "Print the final clique tree";
    "--debug_send", Arg.Unit (fun () -> params.debug_send <- true),
     "Debug: print send_msg information";
    "--time", Arg.Unit (fun () -> params.time <- true),
     "Test time of algorithm";
   ]
  in
  Arg.parse param_specs
    (fun f -> files := !files @ [f]) (* add files *)
    "network_file cpd_file cliquetree_file queries_file";
  match !files with 
  | [nf;cpd;ct;q] -> params.network_file    <- nf;
                     params.cpd_file        <- cpd;
                     params.cliquetree_file <- ct;
                     params.queries_file    <- q
  | _ -> print_endline 
    "Please specify network_file cpd_file cliquetree_file queries_file";
    exit 1

let print_tree tree = print_endline @: string_of_tree tree


let compute_joint cpd_list =
  let {p_of; _} = hd @: parse_queries ~scheme:SumProduct params.queries_file in
  print_endline "applying evidence...";
  let (cpd_list:cpd list) = apply_evidence_to_cpds ~full:false cpd_list p_of in
  print_endline "product...";
  let product = List.fold_left (fun acc cpd ->
      product cpd acc) 
    (empty_cpd ()) cpd_list in
  let cpd = normalize_and_real product in
  let cpd = apply_evidence_to_cpds ~full:true [cpd] p_of in
  print_endline @: string_of_cpd @: hd cpd



let do_inference params cpds query_list =  

  let scheme = if params.action=MaxProductInference 
               then MaxProduct else SumProduct in
  let tree = parse_clique_tree params.cliquetree_file in
  if params.debug_send then print_endline "parsed clique tree";
  set_tree_sepsets tree;
  if params.debug_send then print_endline "set sepsets";
  let cpd_list = cpds in
  if params.debug_send then print_endline "parsed cpds";
  (* debug queries *)
  (*let queries_s = String.concat "\n" @: List.map string_of_query query_list in*)
  (*Printf.printf "queries:\n%s\n" queries_s;*)
  if params.debug_send then print_endline "parsed queries";
  let tree = tree_fill_cpds tree cpd_list in
  save_node_cpds tree;
  if params.debug_send then print_endline "filled tree with cpds";

  let stream_fn tree = 
    upstream tree (fst tree) ~scheme ~print_send:params.debug_send;
    if params.debug_send then print_endline "Downstream...";
    downstream tree ~print_send:params.debug_send;
    if params.print_tree then (print_endline "Tree:"; print_tree tree) else ()
  in
  let answers = 
    if params.debug_send then print_endline "\nWith evidence:";
    process_queries ~incremental:params.incremental stream_fn tree query_list in
  answers

let infer queries ffs num_states num_ts obs : float array = 
  let cpds = Crf.cpds_of_data ffs num_states num_ts obs in
  let p = {
    action=Inference;
    network_file="";
    cpd_file="";
    cliquetree_file="cliquetree.txt";
    queries_file="";
    debug_send=false;
    print_tree=false;
    incremental=true;
    time=false;
  } in
  let answers = List.map unwrap_some @: do_inference p cpds queries in
  Array.of_list answers

