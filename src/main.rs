use serde::{Deserialize, Serialize};
use serde_path_to_error;
use serde_pickle::{
    de::{from_reader, DeOptions},
    Value,
};
use std::io::Read;
use std::{collections::HashMap, fs::File};

use pymol_structs::{
    AtomInfo, Bond, CameraInfoType, ColorRamp, ColorRampElem, ColorRampType, LightInfoType,
    MolInfoType, ObjectInfoType, PSEData, PSEMolecule, SceneInfo, SettingInfo,
};

fn load_pse_data(file_path: &str) -> Result<PSEData, Box<dyn std::error::Error>> {
    let mut file = File::open(file_path)?;
    let mut buffer = Vec::new();
    file.read_to_end(&mut buffer)?;

    let options = DeOptions::new()
        .replace_unresolved_globals()
        .decode_strings();

    let pse_data: PSEData = from_reader(&buffer[..], options)?;
    Ok(pse_data)
}

fn load_and_print_pickle() -> Result<(), Box<dyn std::error::Error>> {
    // Open the file
    // let mut file = File::open("examples/example.pse")?;
    let mut file = File::open("examples/example_molecule_only.pse")?;

    // Read the file contents
    let mut buffer = Vec::new();
    file.read_to_end(&mut buffer)?;

    let options = DeOptions::new()
        .replace_unresolved_globals()
        .decode_strings(); // Decode Python 2 strings as UTF-8
                           // .replace_pyclass_with_dict(true) // Convert Python classes to dictionaries

    // Deserialize the pickle data
    // let deserialized: serde_json::Value = from_reader(&buffer[..], options)?;
    // let deserialized: Value = from_reader(&buffer[..], options)?;
    // println!("Deserialized and re-serialized data:");
    // println!("{}", deserialized);

    let deserialized: PSEData = from_reader(&buffer[..], options)?;
    println!("Deserialized and re-serialized data:");
    println!("{:?}", deserialized);
    // println!("{:?}", deserialized.main[0]);

    // Serialize the data to a pretty-printed JSON string
    // let json_output = serde_json::to_string_pretty(&deserialized)?;

    // Print the serialized output
    // println!("Deserialized and re-serialized data:");
    // println!("{}", json_output);

    // If you want to see the raw pickle serialization, uncomment the following lines:
    // let pickle_output = to_vec(&deserialized)?;
    // println!("\nRaw pickle serialization:");
    // println!("{:?}", pickle_output);
    Ok(())
}

fn main() {
    if let Err(e) = load_and_print_pickle() {
        eprintln!("Error: {}", e);
    }
}
