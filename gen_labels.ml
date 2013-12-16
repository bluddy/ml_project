open Util
open Crf
open Inference

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

let rec sample_until_eof 
  ic ffs num_atoms num_states prev_labels prev_obs num_ts =
  try
   let obs = next_window num_atoms prev_obs ic in
   let labels = sample_next_labels ffs num_states num_ts prev_labels obs in
   let last_label = hd @: List.rev labels in
   print_endline (soi last_label);
   sample_until_eof ic ffs num_atoms num_states labels obs num_ts
  with End_of_file -> ()

let gen_labels file ffs window num_states num_atoms = 
  let ic = open_in file in
  let init_obs = read_n_obs ic window num_atoms in
  let labels = sample_initial_labels ffs num_states window init_obs in
  List.iter (fun state -> print_endline (soi state)) labels;
  sample_until_eof ic ffs num_atoms num_states labels init_obs window;
  close_in ic

let prob_of_lbls obs lbls ffs num_states num_ts = 
  let a = infer ffs num_states num_ts obs in
  (*let p1 = log @: get_p num_ts num_states 1 state a in*)
  let p1 = 0. in (* log 1 *)
  (* look at a pair at a time: prev and cur *)
  let rec prod_rest t = function
    | []  -> failwith "not on empty list"
    | [x] -> 0.
    | l_prev::l_cur::ls ->
      let p = get_p2 num_ts num_states t l_prev l_cur a in
      (log p) +. prod_rest (t+1) (l_cur::ls)
  in exp @: p1 +. (prod_rest 2 @: lbls) 

let next_labels prev_labels lbl_ic = 
        (tl prev_labels) @ [(ios @: input_line lbl_ic)]
