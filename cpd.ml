open Util

(* contains both vars and backptrs *)
type id = int

type cpd_line = id array * float * id array 

type cpd = {vars:id array;
            backptrs: id array;
            data:cpd_line list;
           }

(* save memory with hashtbl *)
let id_str_h = Hashtbl.create 100
let str_id_h = Hashtbl.create 100
let g_cnt = ref 0

let add_id str =
    g_cnt := !g_cnt + 1;
    Hashtbl.add str_id_h str !g_cnt;
    Hashtbl.add id_str_h !g_cnt str;
    !g_cnt

let id_of_str str = try Hashtbl.find str_id_h str 
  with Not_found -> add_id str

let str_of_id id = try Hashtbl.find id_str_h id
 with Not_found -> failwith "Bad id!"

let id_of_str_many strs = List.map id_of_str strs
let str_of_id_many ids = List.map str_of_id ids

let id_of_str_pairs strs =
  List.map (fun (a,b) -> id_of_str a, id_of_str b) strs
let str_of_id_pairs ids =
  List.map (fun (a,b) -> str_of_id a, str_of_id b) ids

let id_of_str_arr strs = Array.map id_of_str strs
let str_of_id_arr ids = Array.map str_of_id ids

let sosa arr = string_of_string_array @: str_of_id_arr arr
let sosl l = string_of_string_list @: str_of_id_many l

let empty_cpd () = {vars=[||]; backptrs=[||]; data=[]}

let cpd_to_log cpd =
  let data = List.rev_map (fun (v,p,back) -> v, log p, back) cpd.data in
  {cpd with data}

let max_cpd_p cpd = match cpd.data with [_,p,_] -> p 
  | _ -> List.fold_left (fun max (_,p,_) -> if p > max then p else max) 
          (let _,p,_ = hd cpd.data in p) 
          (tl cpd.data)

let cpd_from_log cpd =
  match cpd.data with [] -> cpd | _ ->
  let max = max_cpd_p cpd in
  let data = List.rev_map (fun (v,p,b) -> v, exp @: p -. max,b) cpd.data in
  {cpd with data}

let string_of_cpd cpd : string =
  let {vars; backptrs; data} = cpd_from_log cpd in
  let vars = str_of_id_many @: Array.to_list vars in
  let backptrs = str_of_id_many @: Array.to_list backptrs in
    ((match vars with 
     | []  -> "!!! Empty vars !!!"
     | [x] -> x
     | _   -> String.concat ", " vars)^
     (match backptrs with
     | []  -> ""
     | _   -> ": "^String.concat ", " @: backptrs)^
    "\ndata:\n")^
    String.concat "\n" @: 
      List.map (fun (deps, p, back) ->
        let deps_s = string_of_string_array @: str_of_id_arr deps in
        (Printf.sprintf "%s -> %f" deps_s p)^
        match back with [||] -> ""
                      | _    -> ": "^
                         string_of_string_array @: str_of_id_arr back
      ) data

let string_of_cpd_list cs : string = 
  String.concat "\n\n" @: List.map string_of_cpd cs

(* concatenate arrays of strings into one array *)
let concat_vars l_vars =
  let len = List.fold_left (fun acc arr -> acc + Array.length arr)
    0 l_vars in
  let new_arr = Array.create len (-1) in
  let _ = List.fold_left (fun acc_len arr ->
    let arr_len = Array.length arr in
    Array.blit arr 0 new_arr acc_len arr_len;
    acc_len + arr_len
  ) 0 l_vars in
  new_arr

let parse_cpd file =
  (* parse key=val *)
  let get_key_val str = 
    let xs = r_split "=" str in
    hd xs, at xs 1
  in
  let f = read_file file in
  let lines = lines f in
  let h = Hashtbl.create 10 in
  (* accumulate cpds *)
  List.iter (fun line ->
    let elems = r_split " " line in
    let var_name, var_val = get_key_val @: hd elems in
    let var_name, var_val = id_of_str var_name, id_of_str var_val in
    (* get dependencies *)
    let dep_list = r_split "," @: at elems 1 in
    let dep_names, dep_vals = 
      List.fold_right (fun str (acc_names, acc_vals) ->
        let n, v = get_key_val str in
        let n_id, v_id = id_of_str n, id_of_str v in
        n_id::acc_names, v_id::acc_vals) dep_list ([],[]) in
    (* get prob value *)
    let p = float_of_string @: at elems 2 in
    let key = Array.of_list @: var_name::dep_names in
    (* convert to log space already here *)
    let cpd_line = Array.of_list(var_val::dep_vals), log p, [||] in
    let cpd = try Hashtbl.find h key 
      with Not_found -> {vars=key; backptrs=[||]; data=[]} (* create a new cpd *)
    in
    let cpd' = {cpd with data = cpd_line::cpd.data} in
    Hashtbl.replace h key cpd'
  ) lines;
  Hashtbl.fold (fun k v acc -> v::acc) h []

let cpd_find_idxs cpd (var_names:id list) = 
  let h = Hashtbl.create 10 in
  Array.iteri (fun i var -> Hashtbl.add h var i) cpd.vars;
  List.fold_left (fun acc var ->
    Hashtbl.find h var::acc)
    []
    var_names

let cpd_find_idxs_arr cpd (var_names:id array) = 
  let h = Hashtbl.create 10 in
  Array.iteri (fun i var -> Hashtbl.add h var i) cpd.vars;
  Array.fold_left (fun acc var ->
    Hashtbl.find h var::acc)
    []
    var_names

let take_idxs take_idxs len_take xs =
  let take = Array.make len_take (-1) in
  let _ = List.fold_left (fun cnt i ->
    take.(cnt) <- xs.(i); cnt + 1
  ) 0 take_idxs in
  take

(* change from indices taken to indices dropped *)
let invert_idxs idxs count =
  let h = Hashtbl.create 10 in
  List.iter (fun i -> Hashtbl.add h i ()) idxs;
  fst @: iterate (fun (acc,i) ->
      if Hashtbl.mem h i then acc,i+1
      else i::acc, i+1) 
    ([],0) 
    count


(* note: could work straight on indices *)
let marginalize cpd idxs =
  let var_len = Array.length cpd.vars in
  let remain_idxs = invert_idxs idxs var_len in
  let remain_len = List.length remain_idxs in
  let var_names = take_idxs remain_idxs remain_len cpd.vars in
  (* handle data *)
  let h = Hashtbl.create 10 in
  let cpd_real = cpd_from_log cpd in
  List.iter (fun (var_vals, p, _) ->
    let shrunk_vars = take_idxs remain_idxs remain_len var_vals in
    (* check for shrunk vars in hashtable *)
    try 
      (* if we have them, update the probability *)
      let saved_p = Hashtbl.find h shrunk_vars in
      Hashtbl.replace h shrunk_vars @: p +. saved_p
    with Not_found -> (* instantiate *)
      Hashtbl.add h shrunk_vars @: p
  ) cpd_real.data;
  (* get back whole cpd *)
  let data = Hashtbl.fold (fun vars p acc -> (vars, log p, [||])::acc) h [] in
  {cpd with vars=var_names; data}

let marginalize_max cpd idxs =
  let var_len = Array.length cpd.vars in
  let remain_idxs = invert_idxs idxs var_len in
  let remain_len = List.length remain_idxs in
  let var_names = take_idxs remain_idxs remain_len cpd.vars in
  let new_backptrs = take_idxs idxs (var_len-remain_len) cpd.vars in
  let backptrs = concat_vars [cpd.backptrs; new_backptrs] in
  (* handle data *)
  let h = Hashtbl.create 10 in
  List.iter (fun (var_vals, p, old_back_vals) ->
    let shrunk_vars = take_idxs remain_idxs remain_len var_vals in
    (* take care of backpointers *)
    let new_back_vals = take_idxs idxs (var_len-remain_len) var_vals in
    let mixed_back_vals = concat_vars [old_back_vals; new_back_vals] in
    (* check for shrunk vars in hashtable *)
    try 
      (* if we have them, update the probability *)
      let saved_p, _ = Hashtbl.find h shrunk_vars in
      if p > saved_p then
        Hashtbl.replace h shrunk_vars (p, mixed_back_vals);
    with Not_found -> (* instantiate *)
      Hashtbl.add h shrunk_vars (p, mixed_back_vals)
  ) cpd.data;
  (* get back whole cpd *)
  let data = Hashtbl.fold (fun vars (p,back_vals) acc -> 
      (vars, p, back_vals)::acc) 
    h [] in
  {vars=var_names; data; backptrs}

(* filter a cpd by adding evidence, setting a var to a value *)
let add_evidence cpd var_list value_list =
  (*print_endline @: "Adding evidence "^(string_of_string_list var_list)^"\n";*)
  let idxs = cpd_find_idxs cpd var_list in
  let data = List.filter (fun (var_vals, p, _) ->
      List.for_all2 (fun idx v -> var_vals.(idx) = v) idxs value_list)
    cpd.data
  in
  {cpd with data}

(* special intersect returning common indices and diff indices*)
let intersect vars1 vars2 =
  let h = Hashtbl.create 10 in
  Array.iteri (fun i var -> Hashtbl.add h var i) vars1;
  let acc_idx1, acc_idx2, acc_other2 =
    fst @: Array.fold_left 
    (fun ((acc_idx1, acc_idx2, acc_other2),idx2) var2 ->
      try
        let idx1 = Hashtbl.find h var2 in
        Hashtbl.remove h var2;
        (idx1::acc_idx1, idx2::acc_idx2, acc_other2), idx2+1
      with Not_found -> 
        (acc_idx1, acc_idx2, idx2::acc_other2), idx2+1
    ) (([],[],[]),0) vars2
  in
  (* get the diff indices of vars1 using the fact that we removed the
   * common indices *)
  let acc_other1 = Array.fold_left (fun acc_other1 var1 ->
    try
      let idx1 = Hashtbl.find h var1 in
      idx1::acc_other1
    with Not_found -> acc_other1
  ) [] vars1
  in 
  acc_idx1, acc_idx2, acc_other1, acc_other2

(* check that one array is a subset of another *)
let subset vars1 vars2 =
  let _, _, l, _ = intersect vars1 vars2 in
  null l

(* condition a variable to be a certain value and retain the rest of the cpd *)
(* returns a hashtable on the conditioned variables *)
let condition cpd_data idxs_keys idxs_keys_len : (id array, cpd_line list) Hashtbl.t =
  let h = Hashtbl.create 10 in
  match cpd_data with [] -> h | _ ->
  let len_vars = 
    let vars, _, _ = hd cpd_data in
    Array.length vars
  in
  let idxs_vals = invert_idxs idxs_keys len_vars in
  let idxs_vals_len = len_vars - idxs_keys_len in
  List.iter (fun (vars, p, back) ->
    let keys = take_idxs idxs_keys idxs_keys_len vars in
    let vals = take_idxs idxs_vals idxs_vals_len vars in
    try
      let oldvals = Hashtbl.find h keys in
      Hashtbl.replace h keys @: (vals, p, back)::oldvals
    with Not_found -> Hashtbl.add h keys [vals, p, back]
  ) cpd_data;
  h

(* get the join-based product of 2 cpds *)
(* we really use addition since the cpd is in log space *)
(* for backpointers, we don't need to worry about having the same varname because
 * of the running intersection property *)
let product cpd1 cpd2 =
  match cpd1, cpd2 with
  | c, _ when c.data = [] -> cpd2
  | _, c when c.data = [] -> cpd1
  | _, _ ->
  let idxs1, idxs2, diff_idxs1, diff_idxs2  = intersect cpd1.vars cpd2.vars in
  let l_common_idxs = List.length idxs1 in
  let l_diff_idxs1 = List.length diff_idxs1 in
  let l_diff_idxs2 = List.length diff_idxs2 in
  let vars_common = take_idxs idxs1 l_common_idxs cpd1.vars in
  let vars_diff1 = take_idxs diff_idxs1 l_diff_idxs1 cpd1.vars in
  let vars_diff2 = take_idxs diff_idxs2 l_diff_idxs2 cpd2.vars in
  let vars = concat_vars [vars_common; vars_diff1; vars_diff2] in
  let cpd1_c_hash = condition cpd1.data idxs1 l_common_idxs in
  let cpd2_c_hash = condition cpd2.data idxs2 l_common_idxs in
  let backptrs = concat_vars [cpd1.backptrs; cpd2.backptrs] in
  let data =
    Hashtbl.fold (fun keys cpd1_d acc_cpd ->
      let cpd2_d = 
        try Hashtbl.find cpd2_c_hash keys with Not_found -> [] in
      let cross_product =
        List.fold_left (fun acc (vals1, p1, back1) ->
          List.fold_left (fun acc' (vals2, p2, back2) ->
            (concat_vars [keys;vals1;vals2], 
             p1 +. p2, (* log space *)
             concat_vars [back1; back2])::acc'
          ) acc cpd2_d
        ) [] cpd1_d
      in
      List.rev_append cross_product acc_cpd
    ) 
    cpd1_c_hash []
  in
  {vars; backptrs ;data }

(* to do division, we loop over the common vars and divide the bigger factor
 * by the smaller one *)
let div cpd1 cpd2 =
  if cpd2.data = [] then cpd1 else (* handle the '1' for division *)
  let idxs1, idxs2, diff_idxs1, diff_idxs2 = intersect cpd1.vars cpd2.vars in
  if not @: null diff_idxs2 then failwith "Div: second factor too big" else
  let l_common_idxs = List.length idxs1 in
  (* reason to use condition: it rearranges the hash by idx if needed *)
  let cpd2_c_hash = condition cpd2.data idxs2 l_common_idxs in
  (* we mostly keep the first cpd *)
  let new_data =
    List.rev_map (fun (vars, p1, back) ->
        let keys  = take_idxs idxs1 l_common_idxs vars in
        let _, p2, _ = try hd @: Hashtbl.find cpd2_c_hash keys
                    with Not_found -> failwith @: "div: missing value" in
        let div = match p1, p2 with
          (*| 0., 0. -> 0.*)
          (*| _,  0. -> failwith "div: divide by 0"*)
          | _,  _  -> p1 -. p2 (* we're in log space *)
        in
        vars, div, back)
     cpd1.data
  in
  {cpd1 with data=new_data}

let normalize_and_real cpd =
  let cpd_real = cpd_from_log cpd in
  let total_p = List.fold_left (fun acc_p (_,p,_) -> acc_p +. p) 0. cpd_real.data in
  let data' =
    List.fold_left (fun acc (data,p,back) -> 
      (data, p /. total_p, back)::acc) [] cpd_real.data
  in
  {cpd with data=data'}



(* ******* tests *************************************)
(*
let test_cpd = {vars=[|"a";"b"|]; 
                backptrs=[||];
                data=[
                  [|"0";"0"|], 0.5, [||];
                  [|"1";"0"|], 0.5, [||];
                  [|"0";"1"|], 0.2, [||];
                  [|"1";"1"|], 0.8, [||];
                ]
                }

let test_cpd' = {vars=[|"c";"b"|]; 
                backptrs=[||];
                data=[
                  [|"0";"0"|], 0.3, [||];
                  [|"1";"0"|], 0.6, [||];
                  [|"0";"1"|], 0.7, [||];
                  [|"1";"1"|], 0., [||];
                ]
                }

let test_cpd'' = {vars=[|"c";"d"|]; 
                backptrs=[||];
                data=[
                  [|"0";"0"|], 0.3, [||];
                  [|"1";"0"|], 0.6, [||];
                  [|"0";"1"|], 0.7, [||];
                  [|"1";"1"|], 0., [||];
                ]
                }

let condition_test () = condition test_cpd.data [1] 1
  
let product_test () = product test_cpd test_cpd'

let test_cpd2 = {vars=[|"a";"b";"c"|]; 
                backptrs=[||];
                data=[
                  [|"0";"0";"0"|], 0.5, [||];
                  [|"0";"0";"1"|], 0.25, [||];
                  [|"0";"1";"0"|], 0.2, [||];
                  [|"0";"1";"1"|], 0.8, [||];
                  [|"1";"0";"0"|], 0.75, [||];
                  [|"1";"0";"1"|], 0.6, [||];
                  [|"1";"1";"0"|], 0.4, [||];
                  [|"1";"1";"1"|], 0.1, [||];
                ]
                }

let test_cpd2' = {vars=[|"c";"d";"b"|]; 
                backptrs=[||];
                data=[
                  [|"0";"0";"0"|], 0.9, [||];
                  [|"0";"0";"1"|], 0.15, [||];
                  [|"0";"1";"0"|], 0.32, [||];
                  [|"0";"1";"1"|], 0.45, [||];
                  [|"1";"0";"0"|], 0.22, [||];
                  [|"1";"0";"1"|], 0.55, [||];
                  [|"1";"1";"0"|], 0.11, [||];
                  [|"1";"1";"1"|], 0.0, [||];
                ]
                }

  
let condition_test2 () = 
  let h = condition test_cpd2.data [0;2] 2 in
  Hashtbl.find h [|"0";"0"|], 
  Hashtbl.find h [|"0";"1"|], 
  Hashtbl.find h [|"1";"0"|], 
  Hashtbl.find h [|"1";"1"|]
  
let product_test2 () = product test_cpd2 test_cpd2'

let product_test3 () = product test_cpd test_cpd''

let div_test () = div test_cpd2 test_cpd'

let marginalize_test () = marginalize test_cpd2' [0;2]
let marginalize_test2 () = marginalize test_cpd2' [1]

*)
