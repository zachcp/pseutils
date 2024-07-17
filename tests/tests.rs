use pymol_session_utils::molviewspec::nodes::{
    ComponentExpression, ComponentSelector, DownloadParams, KindT, NodeParams, ParseFormatT,
    ParseParams, RepresentationTypeT, State, StructureParams, StructureTypeT,
};
use pymol_session_utils::PSEData;
use serde_json::from_reader;
use std::fs::File;
use std::io::BufReader;
use std::io::Write;

#[test]
fn test_load_pse_data_molecule_only() {
    //https://users.rust-lang.org/t/serde-untagged-enum-ruins-precise-errors/54128/2
    // https://www.gustavwengel.dk/serde-untagged-enum-errors-are-bad
    let deserialized: PSEData = PSEData::load("tests/data/example_molecule_only.pse").unwrap();
    // deserialized.to_json("tests/data/example_molecule_only.json");
    // Check fields of PSEData
    assert!(deserialized.version == 3000000);
}

#[test]
fn test_load_pse_data_molecule_selection() {
    let deserialized: PSEData = PSEData::load("tests/data/example.pse").unwrap();
    assert!(deserialized.version == 3000000);
}

#[test]
fn test_molspecview_json_1cbs() {
    let json_files_component_list = vec![
        "tests/mol-spec-data/1cbs/auth_residue.json",
        "tests/mol-spec-data/1cbs/chain_label.json",
        "tests/mol-spec-data/1cbs/rainbow.json",
        "tests/mol-spec-data/1cbs/validation.json",
    ];

    for json_file in json_files_component_list {
        let file = File::open(json_file).expect(&format!("Failed to open file: {}", json_file));
        let reader = BufReader::new(file);
        let _: Vec<ComponentExpression> = from_reader(reader).expect(&format!(
            "Failed to parse JSON as a Vector of ComponentExpressions: {}",
            json_file
        ));
    }

    // // Todo: Fix the Resideu_range_component/////
    // let json_files_component = vec!["tests/mol-spec-data/1cbs/auth_residue_range.json"];
    // for json_file in json_files_component {
    //     let file = File::open(json_file).expect(&format!("Failed to open file: {}", json_file));
    //     let reader = BufReader::new(file);
    //     let _: ComponentExpression = from_reader(reader).expect(&format!(
    //         "Failed to parse JSON as ComponentExpression: {}",
    //         json_file
    //     ));
    // }
}

#[test]
fn test_molspecview_json_1h9t() {
    let file = File::open("tests/mol-spec-data/1h9t/domains.json").expect("Failed to open file");
    let reader = BufReader::new(file);
    let testvec: Vec<ComponentExpression> =
        from_reader(reader).expect("Failed to parse JSON as ComponentExpression");

    assert_eq!(testvec[0].label_asym_id, Some("A".to_string()));
    assert_eq!(testvec[0].beg_label_seq_id, Some(9));
    assert_eq!(testvec[0].end_label_seq_id, Some(83));

    // todo: these shouw work!
    // assert_eq!(testvec[0].color, "#dd6600");
    // assert_eq!(testvec[0].tooltip, "DNA-binding");

    // let file = File::open("tests/mol-spec-data/1h9t/domains.json").expect("Failed to open file");
    // let reader = BufReader::new(file);
    // let testvec: Vec<ComponentExpression> =
    //     from_reader(reader).expect("Failed to parse JSON as ComponentExpression");
}

#[test]
fn test_molspecview_json_2bvk() {
    let file = File::open("tests/mol-spec-data/2bvk/atoms.json").expect("Failed to open file");
    let reader = BufReader::new(file);
    let testvec: Vec<ComponentExpression> =
        from_reader(reader).expect("Failed to parse JSON as ComponentExpression");

    assert_eq!(testvec[0].atom_index, Some(0));

    // // Todo : Fix color type on Component Expression
    // assert_eq!(testvec[0].color, "#ffdd88");
    // // Todo : Fix tooltip type on Component Expression
    // assert_eq!(testvec[0].tooltip, "First cycle (by atom_index)");
}

// #[test]
// fn test_molspecview_json_full_examples_annotations() {
//     let file =
//         File::open("tests/mol-spec-examples/annotations/state.mvsj").expect("Failed to open file");
//     let reader = BufReader::new(file);
//     let msvj: State = from_reader(reader).expect("Failed to parse JSON as ComponentExpression");

//     // test metadata
//     assert_eq!(msvj.metadata.version, "0.1");
//     assert_eq!(
//         msvj.metadata.title,
//         Some("An example with MVS annotations".to_string())
//     );
//     assert_eq!(
//         msvj.metadata.timestamp,
//         "2024-03-05T18:40:24.870561+00:00".to_string()
//     );

//     // test root and data
//     assert_eq!(msvj.root.kind, KindT::Root);
//     // one child
//     // let download = &msvj.root.children.unwrap()[0];
//     // assert_eq!(download.kind, KindT::Download);
//     // assert_eq!(
//     //     download.params,
//     //     Some(NodeParams::DownloadParams(DownloadParams{"https://files.wwpdb.org/download/1h9t.cif"})),
//     // );
//     //
// }

// #[test]
// fn test_molspecview_json_full_examples_basic() {
//     let file = File::open("tests/mol-spec-examples/basic/state.mvsj").expect("Failed to open file");
//     let reader = BufReader::new(file);
//     let msvj: State = from_reader(reader).expect("Failed to parse JSON as an MVSJ State Object");

//     // test metadata
//     assert_eq!(msvj.metadata.version, "0.1");
//     assert_eq!(
//         msvj.metadata.timestamp,
//         "2023-11-27T12:05:32.145284".to_string()
//     );

//     // test root and data
//     assert_eq!(msvj.root.kind, KindT::Root);
//     // one child
//     let download = &msvj.root.children.unwrap()[0];
//     assert_eq!(download.kind, KindT::Download);
//     // assert_eq!(
//     //     download.params,
//     //     Some(HashMap::from([(
//     //         "url".to_string(),
//     //         serde_json::Value::String("https://files.wwpdb.org/download/1cbs.cif".to_string())
//     //     )])),
//     // );
// }

#[test]
fn test_pdb() {
    let psedata: PSEData = PSEData::load("tests/data/example.pse").unwrap();
    let names = psedata.get_session_names();
    print!("{:?}", names);
    // this has a Molecule and a selection
    assert_eq!(names.len(), 2);

    let mols = psedata.get_molecule_data();
    assert_eq!(mols.len(), 1);
    let atom01 = mols[0].get_atom(0);
    assert!(atom01.x() == 50.87300109863281);
    assert!(atom01.y() == 32.97800064086914);
    assert!(atom01.z() == 2.38700008392334);

    let chains = mols[0].get_chains();
    println!("Chains: {:?}", chains);
    let residues = mols[0].get_residues_by_chain(chains[0].clone());
    println!("Residues: {:?}", residues);

    let residue = mols[0].create_residue(chains[0].clone(), residues[0]);
    println!("Residue: {:?}", residue);

    let chain = mols[0].create_chain(chains[0].clone());
    println!("Chain: {:?}", chain);

    // Check symmetry code
    let (unit, sym) = mols[0].get_unit_cell_symmetry();
    println!("{:?},{:?}", unit, sym);

    // Move on to PDB, baby!
    let pdb = mols[0].to_pdb();

    let _ = pdbtbx::save_pdb(
        &pdb,
        "/Users/zcpowers/Desktop/PSE/pickletest/test_01.pdb",
        pdbtbx::StrictnessLevel::Strict,
    )
    .expect("PDB output");
}

#[test]
fn test_pse_output() {
    // https://colab.research.google.com/drive/1O2TldXlS01s-YgkD9gy87vWsfCBTYuz9#scrollTo=U256gC0Tj2vS
    // builder = mvs.create_builder()
    // (
    //     builder.download(url='https://files.wwpdb.org/download/1cbs.cif')
    //     .parse(format='mmcif')
    //     .assembly_structure(assembly_id='1')
    //     .component()
    //     .representation()
    //     .color(color='#1b9e77')
    // )

    // parse params
    let structfile = ParseParams {
        format: ParseFormatT::Mmcif,
    };

    // struct params
    let structparams = StructureParams {
        structure_type: StructureTypeT::Assembly,
        assembly_id: Some('1'.to_string()),
        ..Default::default()
    };

    // define the component which is going to be `all` here
    let component = ComponentSelector::All("all".to_string());

    // cartoon type
    let cartoon_type = RepresentationTypeT::Cartoon;

    let mut state = State::new();
    state
        .download("https://files.wwpdb.org/download/1cbs.cif")
        .expect("Create a Downlaod node with a URL")
        .parse(structfile)
        .expect("Parseable option")
        .assembly_structure(structparams)
        .expect("a set of Structure options")
        .component(component)
        .expect("defined a valid component")
        .representation(cartoon_type);

    let pretty_json = serde_json::to_string_pretty(&state).unwrap();
    let mut file = File::create("output.json").unwrap();
    file.write_all(pretty_json.as_bytes()).unwrap();
}
