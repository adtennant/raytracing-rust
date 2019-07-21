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
    fn point_at(&self, distance: impl Into<f64>) -> Vector3 {
        self.origin + (self.direction * distance.into())
    }
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
        Sphere { center, radius: radius.into() }
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

    println!("P3");
    println!("{} {}", width, height);
    println!("255");

    let lower_left = Vector3::new(-2.0, -1.0, -1.0);
    let horizontal = Vector3::new(4.0, 0.0, 0.0);
    let vertical = Vector3::new(0.0, 2.0, 0.0);
    let origin = Vector3::new(0.0, 0.0, 0.0);

    let sphere1 = Sphere::new(Vector3::new(0.0, 0.0, -1.0), 0.5);
    let sphere2 = Sphere::new(Vector3::new(0.0, -100.5, -1.0), 100.0);

    let world = World(vec![Box::new(sphere1), Box::new(sphere2)]);

    for y in (0..height).rev() {
        for x in 0..width {
            let u = f64::from(x) / f64::from(width);
            let v = f64::from(y) / f64::from(height);

            let r = Ray {
                origin,
                direction: lower_left + horizontal * u + vertical * v,
            };
            let color = color(&r, &world);

            println!(
                "{} {} {}",
                (color.r * 255.99) as u8,
                (color.g * 255.99) as u8,
                (color.b * 255.99) as u8
            );
        }
    }
}
