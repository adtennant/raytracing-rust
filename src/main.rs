struct Color {
    r: f64,
    g: f64,
    b: f64,
}

impl From<Vec3> for Color {
    fn from(vec: Vec3) -> Self {
        Color {
            r: vec.x,
            g: vec.y,
            b: vec.z,
        }
    }
}

#[derive(Clone)]
struct Vec3 {
    x: f64,
    y: f64,
    z: f64,
}

impl Vec3 {
    fn length(&self) -> f64 {
        (self.x * self.x + self.y * self.y + self.z * self.z).sqrt()
    }

    fn normalized(&self) -> Vec3 {
        self.clone() / self.length()
    }
}

impl std::ops::Add<Vec3> for Vec3 {
    type Output = Self;

    fn add(self, rhs: Vec3) -> Self::Output {
        Vec3 {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
            z: self.z + rhs.z,
        }
    }
}

impl std::ops::Mul<f64> for Vec3 {
    type Output = Self;

    fn mul(self, rhs: f64) -> Self::Output {
        Vec3 {
            x: self.x * rhs,
            y: self.y * rhs,
            z: self.z * rhs,
        }
    }
}

impl std::ops::Div<f64> for Vec3 {
    type Output = Self;

    fn div(self, rhs: f64) -> Self::Output {
        Vec3 {
            x: self.x / rhs,
            y: self.y / rhs,
            z: self.z / rhs,
        }
    }
}

struct Ray {
    origin: Vec3,
    direction: Vec3,
}

fn color(ray: &Ray) -> Color {
    let normalized_direction = ray.direction.normalized();
    let t = 0.5 * (normalized_direction.y + 1.0);

    Color::from(
        Vec3 {
            x: 1.0,
            y: 1.0,
            z: 1.0,
        } * (1.0 - t)
            + Vec3 {
                x: 0.5,
                y: 0.7,
                z: 1.0,
            } * t,
    )
}

fn main() {
    let width = 200;
    let height = 100;

    println!("P3");
    println!("{} {}", width, height);
    println!("255");

    let lower_left = Vec3 {
        x: -2.0,
        y: -1.0,
        z: -1.0,
    };
    let horizontal = Vec3 {
        x: 4.0,
        y: 0.0,
        z: 0.0,
    };
    let vertical = Vec3 {
        x: 0.0,
        y: 2.0,
        z: 0.0,
    };
    let origin = Vec3 {
        x: 0.0,
        y: 0.0,
        z: 0.0,
    };

    for y in (0..height).rev() {
        for x in 0..width {
            let u = f64::from(x) / f64::from(width);
            let v = f64::from(y) / f64::from(height);

            let r = Ray {
                origin: origin.clone(),
                direction: lower_left.clone() + horizontal.clone() * u + vertical.clone() * v,
            };
            let color = color(&r);

            println!(
                "{} {} {}",
                (color.r * 255.99) as u8,
                (color.g * 255.99) as u8,
                (color.b * 255.99) as u8
            );
        }
    }
}
