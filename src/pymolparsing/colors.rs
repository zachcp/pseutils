use once_cell::sync::Lazy;

struct Color {
    name: &'static str, // Changed from String to &'static str,
    r: f32,
    g: f32,
    b: f32,
}

pub static COLOR_SET: Lazy<[Color; 54]> = Lazy::new(|| {
    [
        Color {
            name: "white",
            r: 1.0,
            g: 1.0,
            b: 1.0,
        },
        Color {
            name: "black",
            r: 0.0,
            g: 0.0,
            b: 0.0,
        },
        Color {
            name: "blue",
            r: 0.0,
            g: 0.0,
            b: 1.0,
        },
        Color {
            name: "green",
            r: 0.0,
            g: 1.0,
            b: 0.0,
        },
        Color {
            name: "red",
            r: 1.0,
            g: 0.0,
            b: 0.0,
        },
        Color {
            name: "cyan",
            r: 0.0,
            g: 1.0,
            b: 1.0,
        },
        Color {
            name: "yellow",
            r: 1.0,
            g: 1.0,
            b: 0.0,
        },
        Color {
            name: "dash",
            r: 1.0,
            g: 1.0,
            b: 0.0,
        },
        Color {
            name: "magenta",
            r: 1.0,
            g: 0.0,
            b: 1.0,
        },
        Color {
            name: "salmon",
            r: 1.0,
            g: 0.6,
            b: 0.6,
        },
        Color {
            name: "lime",
            r: 0.5,
            g: 1.0,
            b: 0.5,
        },
        Color {
            name: "slate",
            r: 0.5,
            g: 0.5,
            b: 1.0,
        },
        Color {
            name: "hotpink",
            r: 1.0,
            g: 0.0,
            b: 0.5,
        },
        Color {
            name: "orange",
            r: 1.0,
            g: 0.5,
            b: 0.0,
        },
        Color {
            name: "chartreuse",
            r: 0.5,
            g: 1.0,
            b: 0.0,
        },
        Color {
            name: "limegreen",
            r: 0.0,
            g: 1.0,
            b: 0.5,
        },
        Color {
            name: "purpleblue",
            r: 0.5,
            g: 0.0,
            b: 1.0,
        },
        Color {
            name: "marine",
            r: 0.0,
            g: 0.5,
            b: 1.0,
        },
        Color {
            name: "olive",
            r: 0.77,
            g: 0.7,
            b: 0.0,
        },
        Color {
            name: "purple",
            r: 0.75,
            g: 0.0,
            b: 0.75,
        },
        Color {
            name: "teal",
            r: 0.0,
            g: 0.75,
            b: 0.75,
        },
        Color {
            name: "ruby",
            r: 0.6,
            g: 0.2,
            b: 0.2,
        },
        Color {
            name: "forest",
            r: 0.2,
            g: 0.6,
            b: 0.2,
        },
        Color {
            name: "deepblue",
            r: 0.25,
            g: 0.25,
            b: 0.65,
        },
        Color {
            name: "grey",
            r: 0.5,
            g: 0.5,
            b: 0.5,
        },
        Color {
            name: "gray",
            r: 0.5,
            g: 0.5,
            b: 0.5,
        },
        Color {
            name: "carbon",
            r: 0.2,
            g: 1.0,
            b: 0.2,
        },
        Color {
            name: "nitrogen",
            r: 0.2,
            g: 0.2,
            b: 1.0,
        },
        Color {
            name: "oxygen",
            r: 1.0,
            g: 0.3,
            b: 0.3,
        },
        Color {
            name: "hydrogen",
            r: 0.9,
            g: 0.9,
            b: 0.9,
        },
        Color {
            name: "brightorange",
            r: 1.0,
            g: 0.7,
            b: 0.2,
        },
        Color {
            name: "sulfur",
            r: 0.9,
            g: 0.775,
            b: 0.25,
        },
        Color {
            name: "tv_red",
            r: 1.0,
            g: 0.2,
            b: 0.2,
        },
        Color {
            name: "tv_green",
            r: 0.2,
            g: 1.0,
            b: 0.2,
        },
        Color {
            name: "tv_blue",
            r: 0.3,
            g: 0.3,
            b: 1.0,
        },
        Color {
            name: "tv_yellow",
            r: 1.0,
            g: 1.0,
            b: 0.2,
        },
        Color {
            name: "yelloworange",
            r: 1.0,
            g: 0.87,
            b: 0.37,
        },
        Color {
            name: "tv_orange",
            r: 1.0,
            g: 0.55,
            b: 0.15,
        },
        Color {
            name: "br0",
            r: 0.1,
            g: 0.1,
            b: 1.0,
        },
        Color {
            name: "br1",
            r: 0.2,
            g: 0.1,
            b: 0.9,
        },
        Color {
            name: "br2",
            r: 0.3,
            g: 0.1,
            b: 0.8,
        },
        Color {
            name: "br3",
            r: 0.4,
            g: 0.1,
            b: 0.7,
        },
        Color {
            name: "br4",
            r: 0.5,
            g: 0.1,
            b: 0.6,
        },
        Color {
            name: "br5",
            r: 0.6,
            g: 0.1,
            b: 0.5,
        },
        Color {
            name: "br6",
            r: 0.7,
            g: 0.1,
            b: 0.4,
        },
        Color {
            name: "br7",
            r: 0.8,
            g: 0.1,
            b: 0.3,
        },
        Color {
            name: "br8",
            r: 0.9,
            g: 0.1,
            b: 0.2,
        },
        Color {
            name: "br9",
            r: 1.0,
            g: 0.1,
            b: 0.1,
        },
        Color {
            name: "pink",
            r: 1.0,
            g: 0.65,
            b: 0.85,
        },
        Color {
            name: "firebrick",
            r: 0.698,
            g: 0.13,
            b: 0.13,
        },
        Color {
            name: "chocolate",
            r: 0.555,
            g: 0.222,
            b: 0.111,
        },
        Color {
            name: "brown",
            r: 0.65,
            g: 0.32,
            b: 0.17,
        },
        Color {
            name: "wheat",
            r: 0.99,
            g: 0.82,
            b: 0.65,
        },
        Color {
            name: "violet",
            r: 1.0,
            g: 0.5,
            b: 1.0,
        },
    ]
});

#[test]
fn test_first_color_name() {
    assert_eq!(COLOR_SET[0].name, "white");
    assert_eq!(COLOR_SET[20].name, "teal");
    assert_eq!(COLOR_SET[53].name, "violet");
}
