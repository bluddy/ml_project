open Crf
open Util
open Functions
open Cpd
open Inference
open Gen_labels

type action = GradientAscent | GenLabels | TestInference

type params = {
  mutable input_file : string;
  mutable label_file : string;
  mutable window : int;
  mutable ts : int;
  mutable num_atoms : int;
  mutable num_states : int;
  mutable alpha : float;
  mutable num_iter : int;
  mutable action : action;
  mutable queries : Queries.query list;
  mutable clique_tree : CliqueTree.tree;
  mutable sigma_squared : float option;
  mutable epsilon : float;
  mutable data_type : [`Atoms | `Features];
}

let next_window num_atoms last_obs in_chan : obs array = 
  let next = read_one_obs num_atoms in_chan in
  for i=0 to Array.length last_obs - 2 do
    last_obs.(i) <- last_obs.(i+1)
  done;
  last_obs.(Array.length last_obs - 1) <- next;
  last_obs

(* calculate joint p(ys|x) by multiplying factors over all data *)
let calculate_likelihood_atoms
  obs_file label_file ffs window num_states num_atoms infdata sigma_sq = 
  let obs_ic = open_in obs_file in
  let lbl_ic = open_in label_file in
  let init_obs = read_n_obs obs_ic window num_atoms in 
  let init_lbls = List.map ios @: read_n_lines window lbl_ic in 
  let p1 = prob_of_lbls init_obs init_lbls ffs num_states window infdata in
  (* slide the window and complete calculation *)
  let rec loop acc prev_obs prev_lbls =
    try
      let obs = next_window num_atoms prev_obs obs_ic in
      let lbls = next_labels prev_lbls lbl_ic in
      let p = prob_of_lbls obs lbls ffs num_states window infdata in
      loop (p +. acc) obs lbls
    with End_of_file ->
      close_in obs_ic;
      close_in lbl_ic;
      acc
  in 
  let unreg = loop p1 init_obs init_lbls in
  match sigma_sq with
  | None    -> unreg
  | Some s2 ->
    let sig_factor = 1. /. (2. *. s2) in
    (* regularize the log likelihood *)
    let regularizer = 
      List.fold_left (fun agg ff ->
        agg +. ff.weight *. ff.weight *. sig_factor
      ) 0. ffs
    in unreg -. regularizer

(* calculate joint p(ys|x) by multiplying factors over all data *)
let calculate_likelihood_features
  data labels ffs window num_states infdata sigma_sq = 
  let unreg, _, _ =
    List.fold_left2 (fun (sum_p, win_obs, win_y) obs y ->
      array_shiftl 1 win_obs;
      win_obs.(window - 1) <- obs;
      let win_y' = list_drop 1 @: win_y@[y] in
      let p = prob_of_lbls win_obs win_y' ffs num_states window infdata in
      sum_p +. p, win_obs, win_y'
    ) 
    (* add dummy values since those will be removed *)
    (0., 
     Array.of_list @: hd data::list_take (window-1) data, 
     hd labels::list_take (window-1) labels
    )
    (list_drop (window-1) data)
    (list_drop (window-1) labels)
  in
  match sigma_sq with
  | None    -> unreg
  | Some s2 ->
    let sig_factor = 1. /. (2. *. s2) in
    (* regularize the log likelihood *)
    let regularizer = 
      List.fold_left (fun agg ff ->
        agg +. ff.weight *. ff.weight *. sig_factor
      ) 0. ffs
    in unreg -. regularizer

let gradient1 obs lbls window num_states f fcs ps =
  let lblsa = Array.of_list lbls in 
  let c,e = List.fold_left (fun (curr,exp) t -> 
    let prevs = lblsa.(t-2) in
    let cs    = lblsa.(t-1) in
    let fval = apply_ff f prevs cs obs t in 
    let newcurr = curr +. fval  in
    let p = get_p window num_states t fcs ps in
    let otherf = apply_ff f (-1) fcs obs t in 
    (* debug *)
    (*Printf.eprintf "gradient1: fval(%f), otherf(%f), p(%f)\n" fval otherf p;*)
    let newexp = exp +. p *. otherf in
    (newcurr,newexp)
  ) (0.,0.) (create_range 2 (window-1))
  in
  (* debug *)
  (*Printf.eprintf "gradient1: c(%f), e(%f)\n" c e;*)
  c -. e

let gradient2 obs lbls window num_states f fps fcs ps =
  let lblsa = Array.of_list lbls in 
  let c,e = List.fold_left (fun (curr,exp) t ->
    let prevs = lblsa.(t-2) in
    let cs    = lblsa.(t-1) in
    let fval = apply_ff f prevs cs obs t in 
    let newcurr = curr +. fval  in
    let p = get_p2 window num_states t fps fcs ps in
    let otherf = apply_ff f fps fcs obs t in
    (* debug *)
    (*Printf.eprintf "gradient2: fval(%f), otherf(%f), p(%f)\n" fval otherf p;*)
    let newexp = exp +. p *. otherf in
    (newcurr,newexp)
  ) (0.,0.) (create_range 2 (window-1))
  in
  (* debug *)
  (*Printf.eprintf "gradient2: c(%f), e(%f)\n" c e;*)
  c -. e 

(* return list of gradients for the ffs *)
let gradient_step infdata obs lbls ffs window num_states =
  let probs = infer infdata ffs num_states window obs in
  (* debug *)
  (*Array.iter (fun p -> Printf.eprintf "%f, " p) ps;*)
  (*Printf.eprintf "\n";*)
  let grads = List.map (function
    {fn;prev_state;curr_state;_} ->
      match curr_state, prev_state with
      | None, _            -> failwith "error for now"
      | Some fcs, None     -> 
        gradient1 obs lbls window num_states fn fcs probs
      | Some fcs, Some fps ->
        gradient2 obs lbls window num_states fn fps fcs probs
  ) ffs
  in
  (* debug *)
  (*List.iter (fun p -> Printf.eprintf "%f, " p) grads;*)
  (*Printf.eprintf "\n";*)
  grads

let gradient_sweep_atoms p ffs =
  let obs_ic = open_in p.input_file in
  let lbl_ic = open_in p.label_file in
  let init_obs = read_n_obs obs_ic p.window p.num_atoms in 
  let init_lbls = List.map ios @: read_n_lines p.window lbl_ic in 
  let num_slides = foi @: p.ts - p.window + 1 in
  let inv_num_slides = 1. /. num_slides in
  let grads_init = list_populate (fun _ -> 0.) 0 (List.length ffs)
  in
  let rec loop prev_obs prev_lbls acc_grads =
      let grads = gradient_step (p.queries, p.clique_tree) prev_obs
        prev_lbls ffs p.window p.num_states in
      let acc_grads' = List.map2 (+.) grads acc_grads in
      try
        let obs  = next_window p.num_atoms prev_obs obs_ic in
        let lbls = next_labels prev_lbls lbl_ic in
        loop obs lbls acc_grads'
      with End_of_file ->
        close_in obs_ic;
        close_in lbl_ic;
        acc_grads'
  in
  let grads = loop init_obs init_lbls grads_init in
  let inv_sigma_sq = match p.sigma_squared with
    | None    -> 0.
    | Some s2 -> 1. /. s2
  in
  let ffs' = List.map2 (fun grad ff ->
    let weight = ff.weight +. (grad *. p.alpha *. inv_num_slides) -.
      ff.weight *. inv_sigma_sq in
    (* debug *)
    (*Printf.eprintf "%f, " weight;*)
    {ff with weight}
  ) grads ffs
  in
  (* debug *)
  (*Printf.eprintf "\n";*)
  ffs'

let gradient_sweep_features p ffs data labels =
  let num_slides = foi @: p.ts - p.window + 1 in
  let inv_num_slides = 1. /. num_slides in
  let grads_init = list_populate (fun _ -> 0.) 0 (List.length ffs) in
  let grads, _, _ =
    List.fold_left2 (fun (acc_grads, win_obs, win_y) obs y ->
      array_shiftl 1 win_obs;
      win_obs.(p.window - 1) <- obs;
      let win_y' = list_drop 1 @: win_y@[y] in
      let grads = gradient_step (p.queries, p.clique_tree) win_obs
        win_y' ffs p.window p.num_states in
      let acc_grads' = List.map2 (+.) grads acc_grads in
      acc_grads', win_obs, win_y'
    ) 
    (* add dummy values since those will be removed *)
    (grads_init, 
     Array.of_list @: hd data::list_take (p.window-1) data, 
     hd labels::list_take (p.window-1) labels
    )
    (list_drop (p.window-1) data)
    (list_drop (p.window-1) labels)
  in
  let inv_sigma_sq = match p.sigma_squared with
    | None    -> 0.
    | Some s2 -> 1. /. s2
  in
  let ffs' = List.map2 (fun grad ff ->
    let weight = ff.weight +. (grad *. p.alpha *. inv_num_slides) -.
      ff.weight *. inv_sigma_sq in
    (* debug *)
    (*Printf.eprintf "%f, " weight;*)
    {ff with weight}
  ) grads ffs
  in
  (* debug *)
  (*Printf.eprintf "\n";*)
  ffs'
 
let gradient_ascent p ffs =
  let _, labels, data = match p.data_type with
    | `Atoms -> 0, [], []
    | `Features -> read_data_file p.input_file 1
  in
  let compute_ll ffs = match p.data_type with
    | `Atoms    -> calculate_likelihood_atoms p.input_file
                     p.label_file ffs p.window p.num_states p.num_atoms
                     (p.queries, p.clique_tree) p.sigma_squared
    | `Features -> calculate_likelihood_features data labels ffs p.window
                     p.num_states (p.queries, p.clique_tree)
                     p.sigma_squared
  in
  let init_ll = compute_ll ffs in
  print_endline @: sof init_ll;
  let rec loop ffs last_ll = function
    | 0 -> ffs
    | i ->
      let newffs = match p.data_type with
        | `Atoms    -> gradient_sweep_atoms p ffs
        | `Features -> gradient_sweep_features p ffs data labels
      in
      let ll = compute_ll newffs in
      if abs_float(ll -. last_ll) <= p.epsilon then newffs else begin
      print_endline @: sof ll;
      loop newffs ll (i-1) end
  in loop ffs init_ll p.num_iter 

let sort_ffs ffs =
  List.sort (fun {weight=w1;_} {weight=w2;_} -> 
    if (w1) < (w2) then (-1) else 1
  ) ffs

let params = {
  input_file = "";
  label_file = "";
  window = 5;
  ts = 10;
  num_atoms = 100;
  num_states = 5;
  alpha = 0.001;
  num_iter = 1000;
  action = GradientAscent;
  queries = [];
  clique_tree= CliqueTree.empty_tree ();
  sigma_squared=Some 1.;
  epsilon=10e-16;
  data_type=`Features;
}

let main () =
  let param_specs = Arg.align [
    "--labels", Arg.String (fun l -> params.label_file <- l; params.data_type <- `Atoms),
        "label_file     Set the label file";
    "--window", Arg.Int (fun i -> params.window <- i),
        "timestep window     Set the number of timesteps in a window";
    "--ts", Arg.Int (fun i -> params.ts <- i),
        "timesteps           Set the number of total timesteps";
    "--atoms", Arg.Int (fun i -> params.num_atoms <- i),
        "atoms     Set the number of atoms";
    "--states", Arg.Int (fun i -> params.num_states <- i),
        "states     Set the number of states";
    "--gen_labels", Arg.Unit (fun _ -> params.action <- GenLabels),
        "gen_labels     Generate labels for the data";
    "--test_inference", Arg.Unit (fun _ -> params.action <- TestInference),
        "test inference     Test inference";
    "--alpha", Arg.Float (fun f -> params.alpha <- f),
        "alpha     Set alpha";
    "--iter", Arg.Int (fun i -> params.num_iter <- i),
        "iterations     Set the number of iterations";
    "--sigma", Arg.Float (fun f -> params.sigma_squared <- Some f),
        "sigma square     Set sigma squared";
    "--no-sigma", Arg.Unit (fun _ -> params.sigma_squared <- None),
        "sigma square     Don't use sigma squared";
  ] in
  let usage_msg =
    Printf.sprintf "%s obs_file [options]" Sys.argv.(0) in
  Arg.parse param_specs
    (fun f -> params.input_file <- f)
    usage_msg;

  let p = params in
  if p.input_file = "" then print_endline usage_msg else
  let num_states, window, num_atoms, obs_file, alpha, iter = 
    p.num_states, p.window, p.num_atoms, p.input_file,
    p.alpha, p.num_iter
  in 
  (* generate queries for whoever needs them later *)
  let qs = Queries.gen_queries p.num_states p.window in
  p.queries <- qs;
  (* create a clique tree for inference *)
  let clique_t_s = CliqueTree.clique_tree_string_of_crf p.window in
  p.clique_tree <- CliqueTree.parse_clique_tree clique_t_s;
  CliqueTree.set_tree_sepsets p.clique_tree;
  let s = match params.data_type with 
  | `Features -> "Features mode" | `Atoms -> "Atoms mode" in
  print_endline s;

  match params.action with
  | GradientAscent ->
    (* build feature functions *)
    let ffs =  build_all_fns num_states num_atoms p.data_type in
    let newffs = gradient_ascent params ffs in
    print_ffs @: sort_ffs newffs

  | GenLabels ->
    let ffs = Functions.build_all_fns num_states num_atoms p.data_type in
    gen_labels obs_file ffs window num_states num_atoms (p.queries, p.clique_tree)

  | TestInference ->
    let ffs = build_transition_ffs num_states `Atoms in
    let lambdas=[-0.622;2.59; -2.35; 0.777] in
    let ffs = 
      List.map2 (fun ff weight -> {ff with weight}) ffs lambdas in
    let ps = infer (p.queries, p.clique_tree) ffs 2 3 [||] in
    Array.iter (fun p -> Printf.printf "%f\n" p) ps
    
let _ =
  if !Sys.interactive then ()
  else main ()

