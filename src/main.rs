fn main() {
    let width = 200;
    let height = 100;

    println!("P3");
    println!("{} {}", width, height);
    println!("255");

    for y in (0..height).rev() {
        for x in 0..width {
            let r = f64::from(x) / f64::from(width);
            let g = f64::from(y) / f64::from(height);
            let b = 0.2;

            println!("{} {} {}", (r * 255.99) as u8, (g * 255.99) as u8, (b * 255.99) as u8);
        }
    }
}
