open Util
open Cpd

type state = int
type atom = {
  x : float;
  y : float;
  z : float;
  orig_idx : int;
}

type obs_atom = {
  atoms : atom array;
  timestep : int;
}

type obs_feature = {
  timestep : int;
  chi1 : float array;
  chi2 : float array;
  h_bonds : int array;
}

type fn = 
               (*last,  curr     ,                  , t    , value*)
  | Atom of    (state -> state -> obs_atom array -> int -> float)
  | Feature of (state -> state -> obs_feature array -> int -> float)

type feature_fn = {
  comment: string;
  weight: float;
  atom_idx   : int option;
  prev_state : int option;
  curr_state : int option;
  fn : fn;
}

type obs_features = obs_feature list

let get_p num_ts num_states curr_t curr_s ps =
  let i = (num_states*(curr_t-1)) + (curr_s - 1) in
  ps.(i)

let get_p2 num_ts num_states curr_t prev_s curr_s ps =
  let i = (num_states*num_states*(curr_t-2)) + (num_states*(prev_s - 1)) + (curr_s-1) in
  let num_in_front = num_ts*num_states in
  ps.(num_in_front + i)

let read_one_obs num_atoms ic : obs_atom = 
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
let read_n_obs ic n num_atoms : obs_atom array =
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
      let fn = match fn with
      | Atom f -> f | _ -> failwith "only atom functions for now"
      in
      acc +. (weight *. (fn prev_state curr_state obs t))
    ) 0. ffs 
  (*Parmap.parfold
    ~ncores:2
    (fun {weight;fn;_} acc ->
      acc +. (weight *. (fn prev_state curr_state obs t))
    )
    (Parmap.L ffs)
    0.
    (+.)*)

let random_weight () = (Random.float 1.) -. 0.5

(* feature functions should only be called starting at t=2 *)
let build_1state_xffs num_states num_atoms = []
(*
     List.flatten @:
     list_populate (fun on_state ->
       list_populate (fun atom_num ->
         {
          comment = Printf.sprintf
            "On state %i decreasing atom[%i].x" on_state atom_num;
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
*)

let build_1state_xffs2 num_states num_atoms = []
(*
     List.flatten @:
     list_populate (fun on_state ->
       list_populate (fun atom_num ->
         {
          comment = Printf.sprintf
            "On state %i increasing atom[%i].x" on_state atom_num;
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
*)

let build_transition_ffs num_states = []
(*
  List.flatten @: 
  list_populate (fun prev_on_state ->
    list_populate (fun on_state ->
      {
        comment = Printf.sprintf
            "State transition %i-%i" prev_on_state on_state;
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
*)

let build_1state_cont num_states num_atoms = []
(*
     List.flatten @:
     list_populate (fun on_state ->
       list_populate (fun atom_num ->
         {
          comment = Printf.sprintf
            "On state %i relative change atom[%i].x" on_state atom_num;
          atom_idx = Some atom_num;
          prev_state=None;
          curr_state=Some on_state;
          weight=random_weight ();
          fn = fun last curr obs t ->
                if curr = on_state then
                  let cur_val = obs.(t-2).atoms.(atom_num).x in
                  (obs.(t-1).atoms.(atom_num).x -. cur_val) /. cur_val
                else 0.
         }
       ) 0 num_atoms
     ) 1 num_states 
*)

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
          [|y_last_id; y_curr_id|], p, [||]
        ) 1 num_states
      ) 1 num_states
    in
    let t_last_id = id_of_str @: soi (t-1) in
    let t_id = id_of_str @: soi t in
    {vars=[|t_last_id;t_id|]; backptrs=[||];data}
  ) 2 (num_timesteps-1)

let next_window num_atoms last_obs in_chan : obs_atom array = 
  let next = read_one_obs num_atoms in_chan in
  for i=0 to Array.length last_obs - 2 do
    last_obs.(i) <- last_obs.(i+1)
  done;
  last_obs.(Array.length last_obs - 1) <- next;
  last_obs

(* read all data from a file/string *)
let read_data str begin_ts =
  let r_ts = Str.regexp "{" in
  let r_data = Str.regexp "\\[" in
  let r_comma = Str.regexp ", " in
  let ts = Str.split r_ts str in
  let ts = list_drop 1 ts in (* just a beginning [ *)
  let ts = List.map (fun t -> 
    let data = Str.split r_data t in (* split on [ *)
    let len = List.length data in
    if len <> 5 then failwith @: "Incorrect split length "^soi len else
    let data = list_drop 1 data in
    List.map (fun d ->
      let vals = Str.split r_comma d in
      let vals = list_drop_end 1 vals in (* label for next data *)
      let vals' = List.rev vals in
      let vals' = match vals' with
        | []    -> failwith "Empty list"
        | x::xs -> str_drop_end 1 x::xs (* drop the ] from the number *)
      in
      List.rev vals'
    ) data
  ) ts
  in
  (* convert to arrays of numbers *)
  let ts = List.map (fun t ->
    List.map (fun d -> List.map (fun x -> fos x) d) t) ts
  in
  (* parse into obs *)
  let y_range = create_range 1 5 in
  let next_idx, ys, obs =
    List.fold_left (fun (acc_idx, acc_y, acc_obs) d ->
      let ($) = List.nth in
      let chi2, chi1, h_bonds, rmsd = Array.of_list (d$0), Array.of_list (d$1), d$2, d$3 in
      (* h_bonds should be ints *)
      let h_bonds = Array.of_list @: List.map (fun f -> iof f) h_bonds in
      let rmsd_y = list_zip y_range rmsd in
      let _, y = list_min_op (fst) rmsd_y in
      acc_idx+1, y::acc_y, {timestep=acc_idx;chi1;chi2;h_bonds}::acc_obs
    ) (begin_ts, [], []) ts
  in
  next_idx, List.rev ys, List.rev obs

let read_data_file file begin_ts =
  let f = read_file file in
  read_data f begin_ts

