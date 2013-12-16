open Util
open Cpd

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

let read_one_obs num_atoms ic : observation = 
  let lines = read_n_lines num_atoms ic in
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
  let ts_atom_tup_list  = 
    List.fold_left (fun acc line ->
      let ts, i, x,y,z = read_line line in
      let newatom = {x;y;z;orig_idx=i} in 
      match acc with
      | []  -> [(ts, [newatom])]
      | (old_ts, alist)::xs when ts=old_ts ->
          (ts, newatom::alist)::xs
      | (old_ts,_)::xd ->
          (ts, [newatom])::acc
    ) [] lines
  in match List.length ts_atom_tup_list with
  | 1 ->
      let data = 
        List.map (fun tup ->
          match tup with
          | (ts, atoms) -> {timestep=ts; atoms=Array.of_list (List.rev atoms)}
      ) ts_atom_tup_list
      in (hd data)
 | _ -> failwith "failed observation read"

(* read data into an observation array *)
let read_n_obs ic n num_atoms : observation array =
  let l = list_populate (fun n ->
      read_one_obs num_atoms ic 
  ) 0 n in
  Array.of_list l

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
  List.iter (fun {weight;comment} ->
     Printf.printf "%s: %f\n" comment weight
  ) ffs

let get_potential ffs prev_state curr_state obs t = 
  List.fold_left (fun acc {weight;fn;_}  ->
    acc +. (weight *. (fn prev_state curr_state obs t))
  ) 0. ffs

let random_weight () = (Random.float 1.) -. 0.5

(* feature functions should only be called starting at t=2 *)
let build_1state_xffs num_states num_atoms =
     List.flatten @:
     list_populate (fun on_state ->
       list_populate (fun atom_num ->
         {
          comment = "X coordinate";
          atom_idx = Some atom_num;
          prev_state=None;
          curr_state=Some on_state;
          weight=random_weight ();
          fn = fun last curr obs t ->
                if curr = on_state &&
                   obs.(t-1).atoms.(atom_num).x < obs.(t-2).atoms.(atom_num).x then 1.
                else 0.
         }
       ) 0 num_atoms
     ) 1 num_states 

let build_1state_xffs2 num_states num_atoms =
     List.flatten @:
     list_populate (fun on_state ->
       list_populate (fun atom_num ->
         {
          comment = "X coordinate";
          atom_idx = Some atom_num;
          prev_state=None;
          curr_state=Some on_state;
          weight=random_weight ();
          fn = fun last curr obs t ->
                if curr = on_state &&
                   obs.(t-1).atoms.(atom_num).x >= obs.(t-2).atoms.(atom_num).x then 1.
                else 0.
         }
       ) 0 num_atoms
     ) 1 num_states 

let build_transition_ffs num_states = 
  List.flatten @: 
  list_populate (fun prev_on_state ->
    list_populate (fun on_state ->
      {
        comment = "Transition function";
        atom_idx = None;
        prev_state=Some prev_on_state;
        curr_state=Some on_state;
        weight=random_weight ();
        fn = fun last curr obs t ->
               if last = prev_on_state && curr = on_state 
               then 1.
               else 0.
      }
    ) 1 num_states
  ) 1 num_states 

(* build initial cpds for our data *)
let cpds_of_data ffs num_states num_timesteps obs =   
  list_populate (fun t ->
    let data = List.flatten @:
      list_populate (fun y_last ->
        let y_last_id = id_of_str @: soi y_last in
        (* all (last, cur) combination *)
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

let next_window num_atoms last_obs in_chan : observation array = 
  let next = read_one_obs num_atoms in_chan in
  for i=0 to Array.length last_obs - 2 do
    last_obs.(i) <- last_obs.(i+1)
  done;
  last_obs.(Array.length last_obs - 1) <- next;
  last_obs

