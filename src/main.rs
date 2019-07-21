#[derive(Clone, Copy)]
struct Vector3 {
    x: f64,
    y: f64,
    z: f64,
}

impl Vector3 {
    fn new(x: f64, y: f64, z: f64) -> Self {
        Vector3 { x, y, z }
    }

    fn dot(&self, other: &Vector3) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    fn length(&self) -> f64 {
        (self.x * self.x + self.y * self.y + self.z * self.z).sqrt()
    }

    fn normalized(&self) -> Vector3 {
        *self / self.length()
    }
}

impl std::ops::Add<Vector3> for Vector3 {
    type Output = Self;

    fn add(self, rhs: Vector3) -> Self::Output {
        Vector3 {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
            z: self.z + rhs.z,
        }
    }
}

impl std::ops::Sub<Vector3> for Vector3 {
    type Output = Self;

    fn sub(self, rhs: Vector3) -> Self::Output {
        Vector3 {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
            z: self.z - rhs.z,
        }
    }
}

impl std::ops::Mul<f64> for Vector3 {
    type Output = Self;

    fn mul(self, rhs: f64) -> Self::Output {
        Vector3 {
            x: self.x * rhs,
            y: self.y * rhs,
            z: self.z * rhs,
        }
    }
}

impl std::ops::Div<f64> for Vector3 {
    type Output = Self;

    fn div(self, rhs: f64) -> Self::Output {
        Vector3 {
            x: self.x / rhs,
            y: self.y / rhs,
            z: self.z / rhs,
        }
    }
}

struct Ray {
    origin: Vector3,
    direction: Vector3,
}

fn hit_sphere(ray: &Ray, center: Vector3, radius: f64) -> bool {
    let oc = ray.origin - center;

    let a = ray.direction.dot(&ray.direction);
    let b = 2.0 * oc.dot(&ray.direction);
    let c = oc.dot(&oc) - radius * radius;

    let discriminant = b * b - 4.0 * a * c;
    discriminant > 0.0
}

struct Color {
    r: f64,
    g: f64,
    b: f64,
}

impl From<Vector3> for Color {
    fn from(vec: Vector3) -> Self {
        Color {
            r: vec.x,
            g: vec.y,
            b: vec.z,
        }
    }
}

fn color(ray: &Ray) -> Color {
    if hit_sphere(ray, Vector3::new(0.0, 0.0, -1.0), 0.5) {
        return Color {
            r: 1.0,
            g: 0.0,
            b: 0.0,
        };
    }

    let normalized_direction = ray.direction.normalized();
    let t = 0.5 * (normalized_direction.y + 1.0);

    Color::from(Vector3::new(1.0, 1.0, 1.0) * (1.0 - t) + Vector3::new(0.5, 0.7, 1.0) * t)
}

fn main() {
    let width = 200;
    let height = 100;

    println!("P3");
    println!("{} {}", width, height);
    println!("255");

    let lower_left = Vector3::new(-2.0, -1.0, -1.0);
    let horizontal = Vector3::new(4.0, 0.0, 0.0);
    let vertical = Vector3::new(0.0, 2.0, 0.0);
    let origin = Vector3::new(0.0, 0.0, 0.0);

    for y in (0..height).rev() {
        for x in 0..width {
            let u = f64::from(x) / f64::from(width);
            let v = f64::from(y) / f64::from(height);

            let r = Ray {
                origin,
                direction: lower_left + horizontal * u + vertical * v,
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
