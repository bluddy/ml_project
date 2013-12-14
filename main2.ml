open Util
open Cpd
open Inference

type state = int
type atom = {
  x : float;
  y : float;
  z : float;
  orig_idx : int;
}

type observation = {
  timestep : int;
  atoms : atom array;
}

type feature_fn = {
  comment: string;
  weight: float;
  atom_idx   : int option;
  prev_state : int option;
  curr_state : int option;
       (*last,  curr     ,                  , t    , value*)
  fn : (state -> state -> observation array -> int -> float)
}


type obsersvations = observation list

let get_p num_ts num_states curr_t curr_s ps =
  let i = (num_states*(curr_t-1)) + (curr_s - 1) in
  ps.(i)

let get_p2 num_ts num_states curr_t prev_s curr_s ps =
  let i = (num_states*num_states*(curr_t-2)) + (num_states*(prev_s - 1)) + (curr_s-1) in
  let num_in_front = num_ts*num_states in
  ps.(num_in_front + i)

(* read data into an observation array *)
let read_obs file : observation array =
  let lines = read_file_lines file in
  let r_comma = Str.regexp "," in
  let read_line line = 
    match Str.split r_comma line with
    | [a;b;c;d;e] ->
      let ts,i,x,y,z =
        ios a, ios b, fos c, fos d, fos e
      in
        ts, i, x,y,z
    | _ -> failwith "bad input file"
  in 
  let data = 
    List.fold_left (fun acc line ->
      let ts, i, x,y,z = read_line line in
      let newatom = {x=x;y=y;z=z;orig_idx=i} in 
      match acc with
      | []  -> [(ts, [newatom])]
      | (old_ts, alist)::xs when ts=old_ts ->
          (ts, newatom::alist)::xs
      | (old_ts,_)::xd ->
          (ts, [newatom])::acc
    ) [] lines
  in 
  let data = 
    List.map (fun tup ->
      match tup with
      | (ts, atoms) -> {timestep=ts; atoms=Array.of_list (List.rev atoms)}
    ) data
  in Array.of_list (List.rev data)

let print_obs obs = 
  let data = Array.to_list obs in
  List.iter (fun tup ->
    match tup with
    | {timestep=ts;atoms=aarray} ->
      List.iter (fun atom -> 
        match atom with 
        | {x=a;y=b;z=c;orig_idx=i} ->
            Printf.printf "Timestep %d, Atom %d: %f,%f,%f\n" ts i a b c
      ) (Array.to_list aarray)
  ) data

let print_ffs ffs =
   List.iter (fun ff ->
     match ff with
     |  {weight=w;comment=c} -> Printf.printf "%s: %f\n" c w
   ) ffs


let get_potential ffs prev_state curr_state obs t = 
   List.fold_left (fun acc ff  ->
     match ff with
     | {weight=w;fn=f} -> acc +. (w *. (f prev_state curr_state obs t))
   ) 0. ffs

let random_weight () = (Random.float 1.) -. 0.5

(* feature functions should only be called starting at t=2 *)
let build_1state_xffs num_states num_atoms =
     let nested_list = list_populate (fun on_state ->
       let newfns = list_populate (fun atom_num ->
         {
          comment = "X coordinate";
          atom_idx = Some atom_num;
          prev_state=None;
          curr_state=Some on_state;
          weight=random_weight ();
          fn = (fun last curr obs t ->
                  if curr = on_state then obs.(t-1).atoms.(atom_num).x else 0.)
         }
       ) 0 num_atoms
       in newfns
     ) 1 num_states 
   in List.flatten nested_list  

let build_transition_ffs num_states = 
     let nested_list = list_populate (fun prev_on_state ->
       let newfns = list_populate (fun on_state ->
         {
          comment = "Transition function";
          atom_idx = None;
          prev_state=Some prev_on_state;
          curr_state=Some on_state;
          weight=random_weight ();
          fn = (fun last curr obs t ->
                  if (last = prev_on_state && curr = on_state) 
                  then 1.
                  else 0.)
         }
       ) 1 num_states
       in newfns
     ) 1 num_states 
   in List.flatten nested_list  

let cpds_of_data ffs num_states num_timesteps obs =   
  list_populate (fun t ->
    let data = List.flatten @:
      list_populate (fun y_last ->
        let y_last_id = id_of_str @: soi y_last in
        list_populate (fun y_curr ->
          let y_curr_id = id_of_str @: soi y_curr in
          let p = get_potential ffs y_last y_curr obs t in
          ([|y_last_id;y_curr_id|],p,[||])
        ) 1 num_states
      ) 1 num_states
    in
    let t_last_id = id_of_str @: soi (t-1) in
    let t_id = id_of_str @: soi t in
    {vars=[|t_last_id;t_id|]; backptrs=[||];data}
  ) 2 (num_timesteps-1)  


let to_non_option l = 
  List.map (fun o -> 
    match o with
    | Some f -> f
    | None -> failwith "Error in query result!"
  ) l 


let main () =
  if Array.length Sys.argv <> 2 then
    Printf.printf "%s file\n" Sys.argv.(0)
  else
    let file = Sys.argv.(1) in
    let obs = read_obs file in
    let ffs = (build_1state_xffs 5 5) @ (build_transition_ffs 5) in
    let cpds = cpds_of_data ffs 5 5 obs in
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
    let a = Array.of_list answers in
    List.iter (print_endline) (List.map (sof) answers);
    print_endline @: sof @: get_p2 5 5 5 5 5 a

let _ =
  if !Sys.interactive then ()
  else main ()




