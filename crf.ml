open Util
open Cpd
open Functions

let get_p num_ts num_states curr_t curr_s ps =
  let i = (num_states*(curr_t-1)) + (curr_s - 1) in
  ps.(i)

let get_p2 num_ts num_states curr_t prev_s curr_s ps =
  let i = (num_states*num_states*(curr_t-2)) + (num_states*(prev_s - 1)) + (curr_s-1) in
  let num_in_front = num_ts*num_states in
  ps.(num_in_front + i)

let read_one_obs num_atoms ic = 
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
          | (ts, atoms) ->
              Obs_atom {timestep=ts; atoms=Array.of_list (List.rev atoms)}
      ) ts_atom_tup_list
      in (hd data)
 | _ -> failwith "failed observation read"

(* read data into an observation array *)
let read_n_obs ic n num_atoms =
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
      acc +. (weight *. (apply_ff fn prev_state curr_state obs t))
    ) 0. ffs 

  (*Parmap.parfold
    ~ncores:2
    (fun {weight;fn;_} acc ->
      acc +. (weight *. (fn prev_state curr_state obs t))
    )
    (Parmap.L ffs)
    0.
    (+.)*)


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

let next_window num_atoms last_obs in_chan : (obs array) = 
  let next = read_one_obs num_atoms in_chan in
  for i=0 to Array.length last_obs - 2 do
    last_obs.(i) <- last_obs.(i+1)
  done;
  last_obs.(Array.length last_obs - 1) <- next;
  last_obs

(* convert [-180;180] to [0;360] *)
let convert_degrees neg_deg =
  if neg_deg < 0. then 360. +. neg_deg else neg_deg

(* read all data from a string *)
(* return inverted data *)
let read_data_inv str begin_ts =
  let r_ts = Str.regexp "{" in
  let r_data = Str.regexp "\\[" in
  let r_comma = Str.regexp ", " in
  let r_num = Str.regexp "-?[0-9]+\\(\\.[0-9]*\\)?" in
  let ts = Str.split r_ts str in
  let ts = list_drop 1 ts in (* just a beginning [ *)
  let ts = List.map (fun t -> 
    let data = Str.split r_data t in (* split on [ *)
    let len = List.length data in
    if len <> 5 then failwith @: "Incorrect split length "^soi len else
    let data = list_drop 1 data in
    List.map (fun d ->
      let vals = Str.split r_comma d in
      let vals' = List.rev vals in
      (* don't drop if it's rmsd (last one) *)
      if Str.string_match r_num (hd vals') 0
      then
        let vals' = Str.matched_string (hd vals')::(list_drop 1 vals') in
        List.rev vals'
      else 
        let vals' = list_drop 1 vals' in
        if Str.string_match r_num (hd vals') 0 then
          let vals' = Str.matched_string (hd vals')::(list_drop 1 vals') in
          List.rev vals'
        else
          failwith "couldn't find number"
    ) data
  ) ts
  in
  (* convert to arrays of numbers *)
  let ts = List.map (fun t ->
    List.map (fun d -> List.map (fun x -> fos x) d) t) ts
  in
  (* parse into obs *)
  let y_range = create_range 1 5 in
  List.fold_left (fun (acc_idx, acc_y, acc_obs) d ->
    let ($) = List.nth in
    let conv_all_degs = List.map convert_degrees in
    let chi2, chi1, h_bonds, rmsd =
      Array.of_list @: conv_all_degs (d$0), 
      Array.of_list @: conv_all_degs (d$1), 
      d$2, d$3 in
    (* h_bonds should be ints *)
    let h_bonds = Array.of_list @: List.map (fun f -> iof f) h_bonds in
    let rmsd_y = list_zip y_range rmsd in
    let (y,_), _ = list_min_op (snd) rmsd_y in
    (*print_endline @: soi y; [>debug<]*)
    acc_idx+1, y::acc_y, 
    Obs_feature {timestep=acc_idx;chi1;chi2;h_bonds}::acc_obs
  ) (begin_ts, [], []) ts

(* regular form of reading data *)
let read_data str begin_ts =
  let idx, ys, xs = read_data_inv str begin_ts in
  idx, List.rev ys, List.rev xs

let read_data_file file begin_ts =
  let f = read_file file in
  read_data f begin_ts
  
(* read multiple data files *)
let read_data_files prefix suffix num =
  let r = create_range 0 num in
  let names = List.map (fun i ->
    prefix^soi i^suffix) r
  in
  let idx, ys, obs =
    List.fold_left (fun (idx,acc_y,acc_obs) name ->
      let f = read_file name in
      let idx',ys,obs = read_data_inv f idx in
      idx',ys@acc_y, obs@acc_obs
    ) (0,[],[]) names
  in 
  idx, List.rev ys, List.rev obs

