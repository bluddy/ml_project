open Crf
open Util
open Cpd
open Inference
open Gen_labels

type action = GradientAscent | GenLabels

type params = {
  mutable input_file : string;
  mutable label_file : string;
  mutable window : int;
  mutable ts : int;
  mutable num_atoms : int;
  num_states : int;
  mutable alpha : float;
  mutable num_iter : int;
  mutable action : action;
}

let next_window num_atoms last_obs in_chan : observation array = 
  let next = read_one_obs num_atoms in_chan in
  for i=0 to Array.length last_obs - 2 do
    last_obs.(i) <- last_obs.(i+1)
  done;
  last_obs.(Array.length last_obs - 1) <- next

(* calculate joint p(ys|x) by multiplying factors over all data *)
let calculate_likelihood obs_file label_file ffs window num_states num_atoms = 
  let obs_ic = open_in obs_file in
  let lbl_ic = open_in label_file in
  let init_obs = read_n_obs obs_ic window num_atoms in 
  let init_lbls = List.map ios @: read_n_lines window lbl_ic in 
  let p1 = log @: prob_of_lbls init_obs init_lbls ffs num_states window in
  (* slide the window and complete calculation *)
  let rec loop acc prev_obs prev_lbls =
    try
      let obs = next_window num_atoms prev_obs obs_ic in
      let lbls = next_labels prev_lbls lbl_ic in
      let p = log @: prob_of_lbls obs lbls ffs num_states window in
      loop (p +. acc) obs lbls
    with End_of_file ->
      close_in obs_ic;
      close_in lbl_ic;
      acc
  in 
  loop p1 init_obs init_lbls

let gradient1 obs lbls num_ts num_states f fcs ps =
  let lblsa = Array.of_list lbls in 
  let c,e = List.fold_left (fun (curr,exp) t -> 
      let prevs = lblsa.(t-2) in
      let cs    = lblsa.(t-1) in
      let fval = f prevs cs obs t in 
      let newcurr = curr +. fval  in
      let p = get_p num_ts num_states t fcs ps in
      let newexp = p *. fval
      in (newcurr,newexp)
    ) (0.,0.) (create_range 2 (num_ts-1))
in ((1.+.c) -. (1.+.e)) -. 2.

let gradient2 obs lbls num_ts num_states f fps fcs ps =
  let lblsa = Array.of_list lbls in 
  let c,e = List.fold_left (fun (curr,exp) t ->
      let prevs = lblsa.(t-2) in
      let cs    = lblsa.(t-1) in
      let fval = f prevs cs obs t in 
      let newcurr = curr +. fval  in
      let p = get_p2 num_ts num_states t fps fcs ps in
      let newexp = p *. fval
      in (newcurr,newexp)
    ) (0.,0.) (create_range 2 (num_ts-1))
  in ((1.+.c) -. (1.+.e)) -. 2.

let gradient_step obs lbls ffs num_ts num_atoms num_states num_examples alpha =
  let a = infer ffs num_states num_ts obs in
  let gradient_ff = function
  | {weight=w;fn=f;prev_state=pso;curr_state=cso} as ff ->
    match cso, pso with
    | None, _ -> failwith "error for now"
    | Some fcs, None -> 
      let d = gradient1 obs lbls num_ts num_states f fcs a in
      let new_w = w +. (alpha *. d /. num_examples) in
      { ff with weight=new_w}
    | Some fcs, Some fps ->
      let d = gradient2 obs lbls num_ts num_states f fps fcs a in
      let new_w = w +. (alpha *. d /. num_examples) in
      { ff with weight=new_w}
  in List.map gradient_ff ffs

let gradient_sweep p ffs =
  let obs_ic = open_in p.input_file in
  let lbl_ic = open_in p.label_file in
  let init_obs = read_n_obs obs_ic p.window p.num_atoms in 
  let init_lbls = List.map ios @: read_n_lines p.window lbl_ic in 
  let num_slides = foi @: p.ts - p.window + 1 in
  let rec loop prev_obs prev_lbls prev_ffs  =
    try
      let newffs = gradient_step
        prev_obs prev_lbls prev_ffs p.window p.num_atoms p.num_states num_slides p.alpha in
      let obs = next_window p.num_atoms prev_obs obs_ic in
      let lbls = next_labels prev_lbls lbl_ic in
      loop obs lbls newffs
    with End_of_file ->
      close_in obs_ic;
      close_in lbl_ic;
      prev_ffs
  in
  loop init_obs init_lbls ffs  
 
let gradient_ascent p ffs =
  let rec loop ffs = function
    | 0 -> ()
    | i ->
      let newffs = gradient_sweep p ffs in
      let ll = calculate_likelihood
        p.input_file p.label_file newffs p.window p.num_states p.num_atoms
      in
      print_endline @: sof ll;
      loop newffs (i-1)
  in loop ffs p.num_iter

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
}

let main () =
  let param_specs = Arg.align [
    "--labels", Arg.String (fun l -> params.label_file <- l),
        "label_file     Set the label file";
    "--window", Arg.Int (fun i -> params.window <- i),
        "timestep window     Set the number of timesteps in a window";
    "--ts", Arg.Int (fun i -> params.ts <- i),
        "timesteps           Set the number of total timesteps";
    "--atoms", Arg.Int (fun i -> params.num_atoms <- i),
        "atoms     Set the number of atoms";
    "--gen_labels", Arg.Unit (fun _ -> params.action <- GenLabels),
        "gen_labels     Generate labels for the data";
    "--alpha", Arg.Float (fun f -> params.alpha <- f),
        "alpha     Set alpha";
    "--iter", Arg.Int (fun i -> params.num_iter <- i),
        "iterations     Set the number of iterations";
  ] in
  let usage_msg =
    Printf.sprintf "%s obs_file [options]" Sys.argv.(0) in
  Arg.parse param_specs
    (fun f -> params.input_file <- f)
    usage_msg;

  let p = params in
  if p.input_file = "" then print_endline usage_msg else
  let num_states, window, num_ts, num_atoms, obs_file, alpha, iter = 
    p.num_states, p.window, p.ts, p.num_atoms, p.input_file,
    p.alpha, p.num_iter
  in 
  match params.action with
  | GradientAscent ->
    let label_file = params.label_file in
    if label_file = "" then print_endline usage_msg else
    (* build feature functions *)
    let ffs = build_1state_xffs num_states num_atoms 
            @ build_transition_ffs num_states in

    let ll =
      calculate_likelihood obs_file label_file ffs num_ts num_states num_atoms
    in print_endline @: sof ll;
    gradient_ascent params ffs

  | GenLabels ->
    let ffs = build_1state_xffs num_states num_atoms 
            @ build_transition_ffs num_states
    in
    gen_labels obs_file ffs num_ts num_states num_atoms
    
let _ =
  if !Sys.interactive then ()
  else main ()

