use rand::prelude::*;

#[derive(Clone, Copy)]
struct Vector3 {
    x: f64,
    y: f64,
    z: f64,
}

impl Vector3 {
    fn new(x: impl Into<f64>, y: impl Into<f64>, z: impl Into<f64>) -> Self {
        Vector3 {
            x: x.into(),
            y: y.into(),
            z: z.into(),
        }
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

impl Ray {
    fn new(origin: Vector3, direction: Vector3) -> Self {
        Ray { origin, direction }
    }

    fn point_at(&self, distance: impl Into<f64>) -> Vector3 {
        self.origin + (self.direction * distance.into())
    }
}

struct Color {
    r: f64,
    g: f64,
    b: f64,
}

impl Color {
    fn new(r: impl Into<f64>, g: impl Into<f64>, b: impl Into<f64>) -> Self {
        Color {
            r: r.into(),
            g: g.into(),
            b: b.into(),
        }
    }
}

impl From<Vector3> for Color {
    fn from(vec: Vector3) -> Self {
        Color::new(vec.x, vec.y, vec.z)
    }
}

impl std::ops::Add<Color> for Color {
    type Output = Self;

    fn add(self, rhs: Color) -> Self::Output {
        Color::new(self.r + rhs.r, self.g + rhs.g, self.b + rhs.b)
    }
}

impl std::ops::Div<f64> for Color {
    type Output = Self;

    fn div(self, rhs: f64) -> Self::Output {
        Color::new(self.r / rhs, self.g / rhs, self.b / rhs)
    }
}

struct Hit {
    t: f64,
    p: Vector3,
    normal: Vector3,
}

trait Hitable {
    fn hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<Hit>;
}

struct Sphere {
    center: Vector3,
    radius: f64,
}

impl Sphere {
    fn new(center: Vector3, radius: impl Into<f64>) -> Self {
        Sphere {
            center,
            radius: radius.into(),
        }
    }
}

impl Hitable for Sphere {
    fn hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<Hit> {
        let oc = ray.origin - self.center;

        let a = ray.direction.dot(&ray.direction);
        let b = oc.dot(&ray.direction);
        let c = oc.dot(&oc) - self.radius * self.radius;

        let discriminant = b * b - a * c;

        if discriminant > 0.0 {
            let t = (-b - discriminant.sqrt()) / a;

            if t > t_min && t < t_max {
                let p = ray.point_at(t);

                return Some(Hit {
                    t,
                    p,
                    normal: (p - self.center) / self.radius,
                });
            }

            let t = (-b + discriminant.sqrt()) / a;

            if t > t_min && t < t_max {
                let p = ray.point_at(t);

                return Some(Hit {
                    t,
                    p,
                    normal: (p - self.center) / self.radius,
                });
            }
        }

        None
    }
}

struct World(Vec<Box<Hitable>>);

impl World {
    fn hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<Hit> {
        let mut nearest_t = t_max;
        self.0
            .iter()
            .filter_map(|h| h.hit(ray, t_min, t_max))
            .fold(None, |result, h| {
                if h.t < nearest_t {
                    nearest_t = h.t;
                    Some(h)
                } else {
                    result
                }
            })
    }
}

struct Camera {
    origin: Vector3,
    lower_left: Vector3,
    horizontal: Vector3,
    vertical: Vector3,
}

impl Camera {
    fn new() -> Self {
        Camera {
            origin: Vector3::new(0.0, 0.0, 0.0),
            lower_left: Vector3::new(-2.0, -1.0, -1.0),
            horizontal: Vector3::new(4.0, 0.0, 0.0),
            vertical: Vector3::new(0.0, 2.0, 0.0),
        }
    }

    fn ray(&self, u: f64, v: f64) -> Ray {
        Ray::new(
            self.origin,
            self.lower_left + self.horizontal * u + self.vertical * v - self.origin,
        )
    }
}

fn color(ray: &Ray, world: &World) -> Color {
    if let Some(hit) = world.hit(ray, 0.0, std::f64::MAX) {
        Color::from(Vector3::new(hit.normal.x + 1.0, hit.normal.y + 1.0, hit.normal.z + 1.0) * 0.5)
    } else {
        let normalized_direction = ray.direction.normalized();
        let t = 0.5 * (normalized_direction.y + 1.0);

        Color::from(Vector3::new(1.0, 1.0, 1.0) * (1.0 - t) + Vector3::new(0.5, 0.7, 1.0) * t)
    }
}

fn main() {
    let width = 200;
    let height = 100;
    let num_samples = 100;

    let mut rng = rand::thread_rng();

    println!("P3");
    println!("{} {}", width, height);
    println!("255");

    let sphere1 = Sphere::new(Vector3::new(0.0, 0.0, -1.0), 0.5);
    let sphere2 = Sphere::new(Vector3::new(0.0, -100.5, -1.0), 100.0);

    let world = World(vec![Box::new(sphere1), Box::new(sphere2)]);
    let camera = Camera::new();

    for y in (0..height).rev() {
        for x in 0..width {
            let pixel = (0..num_samples).fold(Color::new(0, 0, 0), |result, _| {
                let u = (f64::from(x) + rng.gen::<f64>()) / f64::from(width);
                let v = (f64::from(y) + rng.gen::<f64>()) / f64::from(height);

                let r = camera.ray(u, v);
                let color = color(&r, &world);

                result + color
            }) / f64::from(num_samples);

            println!(
                "{} {} {}",
                (pixel.r * 255.99) as u8,
                (pixel.g * 255.99) as u8,
                (pixel.b * 255.99) as u8
            );
        }
    }
}
