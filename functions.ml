open Util

type atom = {
  x : float;
  y : float;
  z : float;
  orig_idx : int;
}


type obs_atom = {
  atoms : atom array;
  a_timestep : int;
}

type obs_feature = {
  timestep : int;
  chi1 : float array;
  chi2 : float array;
  h_bonds : int array;
}

type obs = Obs_atom of obs_atom | Obs_feature of obs_feature

let obs_atom = function Obs_atom x -> x | _ -> failwith "not obs_atom"
let obs_feature = function Obs_feature x -> x | _ -> failwith "not obs_feature"

type observations = obs list

type state = int

type fn = 
                  (*last,  curr     ,                  , t    , value*)
  | Atom_fn    of (state -> state -> obs array -> int -> float)
  | Feature_fn of (state -> state -> obs array -> int -> float)

type feature_fn = {
  comment: string;
  weight: float;
  atom_idx   : int option;
  prev_state : int option;
  curr_state : int option;
  fn : fn;
}

let apply_ff ff last_s curr_s obs_arr t = match ff with
  | Atom_fn f    -> f last_s curr_s obs_arr t
  | Feature_fn f -> f last_s curr_s obs_arr t

let random_weight () = (Random.float 1.) -. 0.5

(* feature functions should only be called starting at t=2 *)
let build_1state_xffs num_states num_atoms =
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
          fn = Atom_fn(fun last curr obs t ->
                if curr = on_state &&
                   (obs_atom obs.(t-1)).atoms.(atom_num).x < 
                   (obs_atom obs.(t-2)).atoms.(atom_num).x then 1.
                else 0.)
         }
       ) 0 num_atoms
     ) 1 num_states 

let build_1state_xffs2 num_states num_atoms =
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
          fn = Atom_fn(fun last curr obs t ->
                if curr = on_state &&
                   (obs_atom obs.(t-1)).atoms.(atom_num).x >= 
                   (obs_atom obs.(t-2)).atoms.(atom_num).x then 1.
                else 0.)
         }
       ) 0 num_atoms
     ) 1 num_states 

let build_transition_ffs num_states fn_type =
  List.flatten @: 
  list_populate (fun prev_on_state ->
    list_populate (fun on_state ->
      let fn last curr _ t =
        if last = prev_on_state && curr = on_state 
        then 1. else 0.
      in
      {
        comment = Printf.sprintf
            "State transition %i-%i" prev_on_state on_state;
        atom_idx = None;
        prev_state=Some prev_on_state;
        curr_state=Some on_state;
        weight=random_weight ();
        fn = match fn_type with 
          | `Atoms    -> Atom_fn fn 
          | `Features -> Feature_fn fn
      }
    ) 1 num_states
  ) 1 num_states 

let build_1state_cont num_states num_atoms =
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
          fn = Atom_fn(fun last curr obs t ->
                if curr = on_state then
                  let cur_val = (obs_atom obs.(t-2)).atoms.(atom_num).x in
                  ((obs_atom obs.(t-1)).atoms.(atom_num).x -. cur_val) /. cur_val
                else 0.)
         }
       ) 0 num_atoms
     ) 1 num_states 

(* feature function *)
let build_bin_fn num_states num_angles bins is_chi1 =
  let is_chi_s = match is_chi1 with true -> "chi1" | false -> "chi2" in
     List.flatten @: List.flatten @:
     list_populate (fun on_state ->
       list_populate (fun angle_num ->
         List.rev @: fst @:
         List.fold_left (fun (acc,last_bin) bin ->
          {
            comment = Printf.sprintf
              "On state %i %s[%i] bin[%.1f-%.1f]" 
                on_state is_chi_s angle_num last_bin bin;
            atom_idx=None;
            prev_state=None;
            curr_state=Some on_state;
            weight=random_weight ();
            fn = Feature_fn(fun last curr obs t ->
              if curr = on_state then
                let chi = if is_chi1 then 
                  (obs_feature obs.(t-1)).chi1
                else
                  (obs_feature obs.(t-1)).chi2
                in
                let v = chi.(angle_num) in
                if bin > last_bin then
                  if v >= last_bin && v < bin then
                    1. else 0.
                else
                  if (v >= last_bin && v <= 360.) ||
                    (v >= 0. && v <= bin) then
                      1. else 0.
              else 0.
            )
          }::acc, bin
         )
         ([], hd bins)
         (tl bins)
       ) 0 num_angles
     ) 1 num_states 

let build_hbond_fn num_states num_hbond =
     List.flatten @:
     list_populate (fun on_state ->
       list_populate (fun hbond_num ->
          {
            comment = Printf.sprintf
              "On state %i hbond[%i]" on_state hbond_num;
            atom_idx=None;
            prev_state=None;
            curr_state=Some on_state;
            weight=random_weight ();
            fn = Feature_fn(fun last curr obs t ->
              if curr = on_state then
                if (obs_feature obs.(t-1)).h_bonds.(hbond_num) = 0 then 0.
                else 1.
              else 0.
            )
          }
       ) 0 num_hbond
     ) 1 num_states 

let build_all_fns_atom num_states num_atoms =
      build_1state_xffs num_states num_atoms 
    @ build_1state_xffs2 num_states num_atoms 
    @ build_transition_ffs num_states `Atoms
    (*@ build_1state_cont num_states num_atoms*)

let build_all_fns_feature num_states (num_chi1,num_chi2,num_hbonds) ~use_hbonds =
  build_transition_ffs num_states `Features @
  build_bin_fn num_states num_chi1 [330.; 30.; 90.; 270.; 330.] true @
  build_bin_fn num_states num_chi2 [330.; 30.; 90.; 270.; 330.] false @
  if use_hbonds then build_hbond_fn num_states num_hbonds else []
      


