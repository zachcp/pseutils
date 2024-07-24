use pseutils::pymolparsing::colors::Color;
use pseutils::pymolparsing::parsing::{CustomValue, SettingsEnum};
use pseutils::PSEData;
const TEST_OUTPUT_DIR: &str = "./test_temporary";

#[test]
fn test_load_pse_data_molecule_only() {
    //https://users.rust-lang.org/t/serde-untagged-enum-ruins-precise-errors/54128/2
    // https://www.gustavwengel.dk/serde-untagged-enum-errors-are-bad
    let deserialized: PSEData = PSEData::load("tests/data/example_molecule_only.pse").unwrap();
    // deserialized.to_json("tests/data/example_molecule_only.json");
    assert!(deserialized.version == 3000000);
    assert!(
        deserialized
            .get_setting(SettingsEnum::Orthoscopic)
            .unwrap()
            .value
            == CustomValue::Integer(0)
    );
    // println!("{:?}", deserialized.settings);
}

#[test]
fn test_load_pse_data_molecule_selection() {
    let deserialized: PSEData = PSEData::load("tests/data/example.pse").unwrap();
    // let _ = deserialized.to_json("tests/data/example.json");
    assert!(deserialized.version == 3000000);
}

#[test]
fn test_pdb_00() {
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
        format!("{}/test_01.pdb", TEST_OUTPUT_DIR),
        pdbtbx::StrictnessLevel::Strict,
    )
    .expect("PDB output");
}

#[test]
fn test_pdb_01() {
    let psedata: PSEData = PSEData::load("tests/data/example.pse").unwrap();
    // let _ = psedata.to_disk(TEST_OUTPUT_DIR);
    let _ = psedata.to_disk_full(TEST_OUTPUT_DIR);
    let url = psedata.to_mvsj_url();
    println!("{}", url);
}

#[test]
fn test_colors() {
    // https://github.com/schrodinger/pymol-open-source/blob/master/layer1/Color.cpp#L880
    let spectrum_s: [[f32; 3]; 13] = [
        [1.0, 0.0, 1.0], // magenta - 0
        [0.5, 0.0, 1.0],
        [0.0, 0.0, 1.0], // blue - 166.66
        [0.0, 0.5, 1.0],
        [0.0, 1.0, 1.0], // cyan - 333.33
        [0.0, 1.0, 0.5],
        [0.0, 1.0, 0.0], // green - 500
        [0.5, 1.0, 0.0],
        [1.0, 1.0, 0.0], // yellow - 666.66
        [1.0, 0.5, 0.0],
        [1.0, 0.0, 0.0], // red - 833.33
        [1.0, 0.0, 0.5],
        [1.0, 0.0, 1.0], // magenta - 999
    ];

    const A_DIV: f32 = 83.333333333;

    // Full spectrum (s000-s999)
    for a in 0..1000 {
        let set1 = (a as f32 / A_DIV) as usize;
        let name = format!("s{:03}", a);
        let f = 1.0 - (a as f32 - (set1 as f32 * A_DIV)) / A_DIV;
        let r = f * spectrum_s[set1][0] + (1.0 - f) * spectrum_s[set1 + 1][0];
        let g = f * spectrum_s[set1][1] + (1.0 - f) * spectrum_s[set1 + 1][1];
        let b = f * spectrum_s[set1][2] + (1.0 - f) * spectrum_s[set1 + 1][2];
        println!(
            " Color {{ name: \"{}\", r: {:?}, g: {:?}, b: {:?} }},",
            name, r, g, b
        );
    }

    let spectrum_r: [[f32; 3]; 13] = [
        [1.0, 1.0, 0.0], // yellow - 0
        [0.5, 1.0, 0.0], // chartreuse
        [0.0, 1.0, 0.0], // green - 166.66
        [0.0, 1.0, 0.5], // limegreen
        [0.0, 1.0, 1.0], // cyan - 333.33
        [0.0, 0.5, 1.0], // marine
        [0.0, 0.0, 1.0], // blue - 500
        [0.5, 0.0, 1.0], // purpleblue
        [1.0, 0.0, 1.0], // magenta - 666.66
        [1.0, 0.0, 0.5], // hotpink
        [1.0, 0.0, 0.0], // red - 833.33
        [1.0, 0.5, 0.0], // orange
        [1.0, 1.0, 0.0], // yellow - 999
    ];

    for a in 0..1000 {
        let set1 = (a as f32 / A_DIV) as usize;
        let name = format!("r{:03}", a);
        let f = 1.0 - (a as f32 - (set1 as f32 * A_DIV)) / A_DIV;
        let r = f * spectrum_r[set1][0] + (1.0 - f) * spectrum_r[set1 + 1][0];
        let g = f * spectrum_r[set1][1] + (1.0 - f) * spectrum_r[set1 + 1][1];
        let b = f * spectrum_r[set1][2] + (1.0 - f) * spectrum_r[set1 + 1][2];
        // Assuming reg_named_color is implemented elsewhere
        println!(
            " Color {{ name: \"{}\", r: {:?}, g: {:?}, b: {:?} }},",
            name, r, g, b
        );
    }

    let spectrum_c: [[f32; 3]; 13] = [
        [1.0, 1.0, 0.0], // yellow - 0
        [0.0, 0.0, 1.0], // blue - 83.333
        [1.0, 0.0, 0.0], // red - 167.67
        [0.0, 1.0, 0.0], // green - 250.00
        [1.0, 0.0, 1.0], // magenta - 333.33
        [0.0, 1.0, 1.0], // cyan - 416.67
        [1.0, 1.0, 0.0], // yellow - 500.00
        [0.0, 1.0, 0.0], // green - 583.33
        [0.0, 0.0, 1.0], // blue - 666.67
        [1.0, 0.0, 1.0], // magenta - 750.00
        [1.0, 1.0, 0.0], // yellow - 833.33
        [1.0, 0.0, 0.0], // red - 916.67
        [0.0, 1.0, 1.0], // cyan - 999
    ];

    let mut name = String::from("c000");
    for a in 0..1000 {
        let set1 = (a as f32 / A_DIV) as usize;
        // sprintf(color->Name,"c%03d",a);
        name = format!("c{:03}", a);
        let f = 1.0 - (a as f32 - (set1 as f32 * A_DIV)) / A_DIV;
        let r = f * spectrum_c[set1][0] + (1.0 - f) * spectrum_c[set1 + 1][0];
        let g = f * spectrum_c[set1][1] + (1.0 - f) * spectrum_c[set1 + 1][1];
        let b = f * spectrum_c[set1][2] + (1.0 - f) * spectrum_c[set1 + 1][2];
        println!(
            " Color {{ name: \"{}\", r: {:?}, g: {:?}, b: {:?} }},",
            name, r, g, b
        );
    }
}
