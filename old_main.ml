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


let compute_joint () =
  let cpd_list = parse_cpd params.cpd_file in
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


let run () = 
  if params.action = ComputeJoint then compute_joint () else ();
  let scheme = if params.action=MaxProductInference 
               then MaxProduct else SumProduct in
  let tree = parse_clique_tree params.cliquetree_file in
  print_endline "parsed clique tree";
  set_tree_sepsets tree;
  print_endline "set sepsets";
  let cpd_list = parse_cpd params.cpd_file in
  print_endline "parsed cpds";
  let query_list = parse_queries ~scheme params.queries_file in
  print_endline "parsed queries";
  let tree = tree_fill_cpds tree cpd_list in
  save_node_cpds tree;
  print_endline "filled tree with cpds";
  match params.action with
  | Print_CPDs -> print_endline @: string_of_cpd_list cpd_list
  | Print_Tree -> print_tree tree
  | Inference  -> 
      let p_time1 = Unix.times () in
      let stream_fn tree = 
        upstream tree (fst tree) ~scheme ~print_send:params.debug_send;
        if params.debug_send then print_endline "Downstream...";
        downstream tree ~print_send:params.debug_send;
        if params.print_tree then (print_endline "Tree:"; print_tree tree) else ()
      in
      stream_fn tree;
      let answers = 
        if params.debug_send then print_endline "\nWith evidence:";
        process_queries ~incremental:params.incremental stream_fn tree query_list in
      let p_time2 = Unix.times () in
      List.iter (function Some a -> Printf.printf "%.13f\n" a
                         | None  -> print_endline "error") answers;
      if params.time then 
        Printf.printf "User time: %f\n" (p_time2.tms_utime -. p_time1.tms_utime)
        else ()

  | MaxProductInference ->
      let stream_fn tree root = 
        upstream tree root ~scheme:MaxProduct ~print_send:params.debug_send;
        if params.print_tree then (print_endline "Tree:"; print_tree tree) else ()
      in
      (* don't do an early pass *)
      let answers = process_queries_max stream_fn tree query_list in
      List.iter (function (p, var_list) ->
          let var_list = str_of_id_pairs var_list in
          Printf.printf "%.13f: %s\n" p (string_of_assignments var_list)
        )
        answers
   | _ -> failwith "bad input"

let _ =
  if !Sys.interactive then ()
  else
    (parse_cmd_line ();
    run ();)


