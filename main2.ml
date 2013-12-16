open Crf
open Util
open Cpd
open Inference


let infer ffs num_states num_ts obs : float array = 
  let cpds = cpds_of_data ffs num_states num_ts obs in
  let p = {
    action=Inference;
    network_file="";
    cpd_file="";
    cliquetree_file="cliquetree.txt";
    queries_file="queries.txt";
    debug_send=false;
    print_tree=false;
    incremental=true;
    time=false;
  } in
  let answers = (to_non_option @: do_inference p cpds) in
  Array.of_list answers

let next_window num_atoms last_obs in_chan : observation array = 
  let next = read_one_obs num_atoms in_chan in
  let l = Array.to_list last_obs in
  let new_l = (tl l) @ [next] in
  Array.of_list new_l

let rec which_state r cdf state : state = 
  if r <= (hd cdf)
  then state
  else which_state r (tl cdf) (state + 1) 



let sample_state ps : state = 
   let r = Random.float 1. in
   let cdf = List.rev @: List.fold_left (fun acc p -> 
      match acc with
      | []    -> [p]
      | (x::xs) -> (p +. x)::x::xs
    ) [] ps in
    which_state r cdf 1

let normalize_list l = 
  let sum = List.fold_left (fun acc p -> p +. acc) 0. l in
  let l = List.map (fun p -> p /. sum) l in
  l

let print_labels labels =
  let rec p n labels = match labels with
  | [] -> ()
  | x::xs -> Printf.printf "Random label for t%d: %d\n" n x;
             p (n+1) (tl labels)
  in p 1 labels
  
let sample_next_state num_states num_ts curr_t last_state a = 
  let pdf = normalize_list @: list_populate (fun n ->
      get_p2 num_ts num_states curr_t last_state n a
  ) 1 num_states in
  sample_state pdf

 
let sample_initial_labels ffs num_states num_ts init_obs = 
  let a = infer ffs num_states num_ts init_obs in
  (* s1 *)
  let s1pdf = normalize_list @: list_populate (fun n ->
      get_p num_ts num_states 1 n a
  ) 1 num_states in
  let s1 = sample_state s1pdf in
  List.rev @: List.fold_left (fun acc t -> match acc with
    | [] -> failwith "put s1 in the acc!"
    | prev_state::xs -> (sample_next_state num_states num_ts t prev_state a)::acc
  ) [s1] (create_range 2 (num_ts -1))
  
let sample_next_labels ffs num_states num_ts prev_labels new_obs =
  let a = infer ffs num_states num_ts new_obs in
  let labels = List.rev @: tl prev_labels in
  let newlast = sample_next_state num_states num_ts num_ts (hd labels) a in
  List.rev @: newlast::labels

let rec sample_until_eof ic ffs num_atoms num_states prev_labels prev_obs num_ts =
  try
   let obs = next_window num_atoms prev_obs ic in
   let labels = sample_next_labels ffs num_states num_ts prev_labels obs in
   let last_label = hd @: List.rev labels in
   print_endline (soi last_label);
   sample_until_eof ic ffs num_atoms num_states labels obs num_ts
  with End_of_file ->
    ()
let gen_labels file ffs num_ts num_states num_atoms = 
    let ic = open_in file in
    let init_obs = read_n_obs ic num_ts num_atoms in
    let labels = sample_initial_labels ffs num_states num_ts init_obs in
    List.iter (fun state -> print_endline (soi state)) labels;
    sample_until_eof ic ffs num_atoms num_states labels init_obs num_ts;
    close_in ic

let prob_of_lbls obs lbls ffs num_states num_ts = 
  let a = infer ffs num_states num_ts obs in
  let state =  hd lbls in
  let p1 = log @: get_p num_ts num_states 1 state a in
  let rec prod_rest t lbls = match lbls with
    | [] -> failwith "not on empty list"
    | [x] -> 0.
    | prev_lbl::lbls ->
      let lbls = tl lbls in
      let curr_lbl = hd lbls in
      let p = get_p2 num_ts num_states t prev_lbl curr_lbl a in
      (log p) +. prod_rest (t+1) lbls
  in exp @: p1 +. (prod_rest 2 lbls) 

let next_labels prev_labels lbl_ic = 
        (tl prev_labels) @ [(ios @: input_line lbl_ic)]


let calculate_likelihood obs_file label_file ffs num_ts num_states num_atoms = 
  let obs_ic = open_in obs_file in
  let lbl_ic = open_in label_file in
  let init_obs = read_n_obs obs_ic num_ts num_atoms in 
  let init_lbls = List.map ios @: read_n_lines num_ts lbl_ic in 
  let p1 = log @: prob_of_lbls init_obs init_lbls ffs num_states num_ts in
  let rec calc_all prev_obs prev_lbls =
    try
      let obs = next_window num_atoms prev_obs obs_ic in
      let lbls = next_labels prev_lbls lbl_ic in
      let p = log @: prob_of_lbls obs lbls ffs num_states num_ts in
      p +. (calc_all obs lbls)
    with End_of_file ->
      close_in obs_ic;
      close_in lbl_ic;
      p1
  in 
  let p = calc_all init_obs init_lbls in
  print_endline (sof p)


let gradient1 obs lbls num_ts num_states f fcs ps =
  let lblsa = Array.of_list lbls in 
  let (c,e) = List.fold_left (fun acc t -> match acc with
     | (curr,exp) ->
       let prevs = lblsa.(t-2) in
       let cs    = lblsa.(t-1) in
       let fval = (f prevs cs obs t) in 
       let newcurr = curr +. fval  in
       let p = (get_p num_ts num_states t fcs ps) in
       let newexp = p *. fval
       in (newcurr,newexp)
     ) (0.,0.) (create_range 2 (num_ts-1))
in ((1.+.c) -. (1.+.e)) -. 2.




let gradient2 obs lbls num_ts num_states f fps fcs ps =
  let lblsa = Array.of_list lbls in 
  let (c,e) = List.fold_left (fun acc t -> match acc with
     | (curr,exp) ->
       let prevs = lblsa.(t-2) in
       let cs    = lblsa.(t-1) in
       let fval = (f prevs cs obs t) in 
       let newcurr = curr +. fval  in
       let p = (get_p2 num_ts num_states t fps fcs ps) in
       let newexp = p *. fval
       in (newcurr,newexp)
     ) (0.,0.) (create_range 2 (num_ts-1))
  in ((1.+.c) -. (1.+.e)) -. 2.


let gradient_step obs lbls ffs num_ts num_atoms num_states num_examples alpha =
  let a = infer ffs num_states num_ts obs in
  let gradient_ff  ff = match ff with
  | {weight=w;fn=f;prev_state=pso;curr_state=cso} ->
    match cso with
    | None -> failwith "error for now"
    | Some fcs ->
      match pso with
      | None -> 
          let d = gradient1 obs lbls num_ts num_states f fcs a in
          let pen = w /. (10. *. (foi num_ts)) in
          let d = d -. pen in 
          let new_w = w +. (alpha *. d /. num_examples) in
          { ff with weight=new_w}
      | Some fps ->
          let d= gradient2 obs lbls num_ts num_states f fps fcs a in
          let pen = w /. (10. *. (foi num_ts)) in
          let d = d -. pen in 
          let new_w = w +. (alpha *. d /. num_examples) in
          { ff with weight=new_w}
  in List.map gradient_ff  ffs

let gradient_sweep obs_file label_file ffs num_ts num_states num_atoms
num_examples alpha = 
  let obs_ic = open_in obs_file in
  let lbl_ic = open_in label_file in
  let init_obs = read_n_obs obs_ic num_ts num_atoms in 
  let init_lbls = List.map ios @: read_n_lines num_ts lbl_ic in 
  let rec calc_all prev_obs prev_lbls prev_ffs  =
    try
      let newffs = gradient_step prev_obs prev_lbls prev_ffs num_ts num_atoms num_states
       num_examples alpha in
      let obs = next_window num_atoms prev_obs obs_ic in
      let lbls = next_labels prev_lbls lbl_ic in
      calc_all obs lbls newffs
    with End_of_file ->
      close_in obs_ic;
      close_in lbl_ic;
      prev_ffs
  in  calc_all init_obs init_lbls ffs  
 
let rec gradient_ascent obs_file label_file ffs num_ts num_states num_atoms niter
nexamples alpha = 
  match niter with
  | 0 -> ()
  | _ ->
    let newffs = gradient_sweep obs_file label_file ffs num_ts num_states
    num_atoms nexamples alpha
    in
    calculate_likelihood obs_file label_file newffs num_ts num_states num_atoms;
    gradient_ascent obs_file label_file newffs num_ts num_states num_atoms
    (niter-1) nexamples alpha

let main () =
  if Array.length Sys.argv <> 3 then
    Printf.printf "%s obs_file label_file\n" Sys.argv.(0)
  else
    let num_ts = 5 in 
    let num_states = 5 in
    let num_atoms = 100 in 
    let obs_file = Sys.argv.(1) in
    let label_file = Sys.argv.(2) in
    let ffs = (build_1state_xffs num_states num_atoms) 
            @ (build_transition_ffs num_states) in
  
    calculate_likelihood obs_file label_file ffs num_ts num_states num_atoms;
    gradient_ascent obs_file label_file ffs num_ts num_states num_atoms 1000
    100. 0.001 
    (*gen_labels obs_file ffs num_ts num_states num_atoms*)
let _ =
  if !Sys.interactive then ()
  else main ()

