struct Color {
    r: f64.
    g: f64,
    b: f64.
}

fn main() {
    let width = 200;
    let height = 100;

    println!("P3");
    println!("{} {}", width, height);
    println!("255");

    for y in (0..height).rev() {
        for x in 0..width {
            let color = Color {
                r: f64::from(x) / f64::from(width),
                g: f64::from(y) / f64::from(height),
                b: 0.2,
            };

            println!(
                "{} {} {}",
                (color.r * 255.99) as u8,
                (color.g * 255.99) as u8,
                (color.b * 255.99) as u8
            );
        }
    }
}
