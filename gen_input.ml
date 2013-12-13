open Util

let gen_data num_atoms num_timesteps =
  List.flatten @:
  list_populate (fun t ->
    list_populate
      (fun i -> t, i, Random.float 1., Random.float 1., Random.float 1.)
      0 num_atoms
  ) 1 num_timesteps

let print_data data =
  List.iter (fun (t, atom, x, y, z) -> 
    Printf.printf "%d,%d,%f,%f,%f\n" t atom x y z
  ) data

let main () =
  if Array.length Sys.argv <> 3 then
    Printf.printf "%s num_atoms num_timesteps\n" Sys.argv.(0)
  else
    let num_atoms = ios @: Sys.argv.(1) in
    let num_timesteps = ios @: Sys.argv.(2) in
    let data = gen_data num_atoms num_timesteps in
    print_data data

let _ = main ()

