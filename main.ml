open Util
open Cpd

type atom = {
  loc : (float * float * float);
}

type observation = {
  timestep : int;
  atoms : atom array;
}

let y_max = 5 (* max state *)

type observations = observation array

type state = int
              (* current, last state                  t     factor *)
type feature_fn = (state -> state -> observations -> int -> float)

(* for now, assume we can read everything into memory *)
let read_data file : observation array =
  let lines = read_file_lines file in
  let r_comma = Str.regexp "," in
  let read_line line = 
    match Str.split r_comma line with
    | [ts;idx;x;y;z] ->
        let timestep, idx, x, y, z =
          ios ts, ios idx, fos x, fos y, fos z
        in
        timestep, idx, x,y,z
    | _ -> failwith "bad input file"
  in
  let data =
    List.fold_left (fun acc line ->
    let ts, idx, x,y,z = read_line line in
    match acc with
      | []              -> [ts, [idx, x,y,z]]
      | (old_ts, atoms)::xs when ts=old_ts -> 
          (ts, (idx, x,y,z)::atoms)::xs
      | (old_ts, _)::xs -> (ts, [idx, x,y,z])::acc
    ) [] lines
  in
  let data = List.rev_map (fun (ts, atoms) ->
    (* sanitize indices and reverse atom list *)
    let _, atoms = List.fold_left (fun (acc_idx, acc) (idx,x,y,z) ->
      (*if acc_idx <> idx then failwith @: "Bad index "^soi idx else*)
        (acc_idx + 1, {loc=(x,y,z)}::acc)) (0,[]) atoms 
    in
    let atoms = Array.of_list atoms in
    {timestep=ts; atoms}
  ) data
  in
  Array.of_list data

(* cur_t: look at current x or the last x *)
(* note: only run from second observation onward *)
let basic_fn cur_t field_num current_y atom_num = 
  fun y_val y_last observations t ->
    let obs = if cur_t then observations.(t) else observations.(t-1) in
    let x, y, z = obs.atoms.(atom_num).loc in
    let v = match field_num with
      | 1 -> x | 2 -> y | 3 -> z
      | _ -> invalid_arg @: "Can't handle "^soi field_num
    in
    if y_val = current_y then v else 0.

let basic_fn_list () : feature_fn list = 
  let t_f = [true; false] in
  let field_num = [1; 2; 3] in
  let current_y = [1; 2] in
  let atom_num = create_range 0 10 in
  let ( * ) = cartesian_product in
  let params = t_f * field_num * current_y * atom_num in
  List.map (fun (((a,b),c),d) -> basic_fn a b c d) params

let cpds_of_data fn_list weight_list (start_ts, end_ts) obs =
  let ts_num = end_ts - start_ts + 1 in
  list_populate (fun t ->
    let data = List.flatten @:
      list_populate (fun y_cur ->
        let y_cur_id = id_of_str @: soi y_cur in
        list_populate (fun y_last ->
          (* calculate factor *)
          let factor = List.fold_left2 (fun acc w fn ->
            acc +. w *. fn y_cur y_last obs t
          ) 0. weight_list fn_list
          in
          let y_last_id = id_of_str @: soi y_last in
          ([|y_last_id; y_cur_id|], factor, [||])
        ) 0 y_max
      ) 0 y_max
    in
    let t_last_id = id_of_str @: soi (t-1) in
    let t_id = id_of_str @: soi t in
    {vars=[|t_last_id; t_id|]; backptrs=[||]; data}
  ) start_ts ts_num

let random_weights len =
  list_populate (fun _ -> (Random.float 1.) -. 0.5) 0 len

let main () =
  if Array.length Sys.argv <> 2 then
    Printf.printf "%s file\n" Sys.argv.(0)
  else
    let file = Sys.argv.(1) in
    let data = read_data file in
    let fn_list = basic_fn_list () in
    let weight_list = random_weights (List.length fn_list) in
    let cpds = cpds_of_data fn_list weight_list (1, 5) data in

let _ =
  if !Sys.interactive then ()
  else main ()
