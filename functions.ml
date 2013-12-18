open Util

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

type obs_features = obs_feature list


type state = int

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
          fn = Atom(fun last curr obs t ->
                if curr = on_state &&
                   obs.(t-1).atoms.(atom_num).x < obs.(t-2).atoms.(atom_num).x then 1.
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
          fn = Atom(fun last curr obs t ->
                if curr = on_state &&
                   obs.(t-1).atoms.(atom_num).x >= obs.(t-2).atoms.(atom_num).x then 1.
                else 0.)
         }
       ) 0 num_atoms
     ) 1 num_states 

let build_transition_ffs num_states =
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
        fn = Atom(fun last curr obs t ->
               if last = prev_on_state && curr = on_state 
               then 1.
               else 0.)
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
          fn = Atom(fun last curr obs t ->
                if curr = on_state then
                  let cur_val = obs.(t-2).atoms.(atom_num).x in
                  (obs.(t-1).atoms.(atom_num).x -. cur_val) /. cur_val
                else 0.)
         }
       ) 0 num_atoms
     ) 1 num_states 

let build_all_fns num_states num_atoms = function
  | `Atoms ->
      build_1state_xffs num_states num_atoms 
    @ build_1state_xffs2 num_states num_atoms 
    @ build_transition_ffs num_states
    (*@ build_1state_cont num_states num_atoms*)
  | `Features -> []

