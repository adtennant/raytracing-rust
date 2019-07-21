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
        self.squared_length().sqrt()
    }

    fn normalized(&self) -> Vector3 {
        *self / self.length()
    }

    fn reflect(&self, n: &Vector3) -> Vector3 {
        *self - *n * self.dot(n) * 2.0
    }

    fn squared_length(&self) -> f64 {
        self.x * self.x + self.y * self.y + self.z * self.z
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

impl std::ops::Mul<Vector3> for Vector3 {
    type Output = Self;

    fn mul(self, rhs: Vector3) -> Self::Output {
        Vector3 {
            x: self.x * rhs.x,
            y: self.y * rhs.y,
            z: self.z * rhs.z,
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

impl std::ops::Mul<f64> for Color {
    type Output = Self;

    fn mul(self, rhs: f64) -> Self::Output {
        Color::new(self.r * rhs, self.g * rhs, self.b * rhs)
    }
}

impl std::ops::Div<f64> for Color {
    type Output = Self;

    fn div(self, rhs: f64) -> Self::Output {
        Color::new(self.r / rhs, self.g / rhs, self.b / rhs)
    }
}

#[derive(Clone, Copy)]
enum Material {
    Lambertian { albedo: Vector3 },
    Metal { albedo: Vector3, fuzz: f64 },
}

impl Material {
    fn scatter(&self, ray: &Ray, hit: &Hit) -> Option<(Ray, Vector3)> {
        match self {
            Material::Lambertian { albedo } => {
                let target = hit.p + hit.normal + random_in_unit_sphere();
                let scattered = Ray::new(hit.p, target - hit.p);
                Some((scattered, *albedo))
            }
            Material::Metal { albedo, fuzz } => {
                let reflected = ray.direction.normalized().reflect(&hit.normal);
                let scattered = Ray::new(hit.p, reflected + random_in_unit_sphere() * *fuzz);

                if scattered.direction.dot(&hit.normal) > 0.0 {
                    Some((scattered, *albedo))
                } else {
                    None
                }
            }
        }
    }
}

struct Hit {
    t: f64,
    p: Vector3,
    normal: Vector3,
    material: Material,
}

trait Hitable {
    fn hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<Hit>;
}

struct Sphere {
    center: Vector3,
    radius: f64,
    material: Material,
}

impl Sphere {
    fn new(center: Vector3, radius: impl Into<f64>, material: Material) -> Self {
        Sphere {
            center,
            radius: radius.into(),
            material,
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
                    material: self.material,
                });
            }

            let t = (-b + discriminant.sqrt()) / a;

            if t > t_min && t < t_max {
                let p = ray.point_at(t);

                return Some(Hit {
                    t,
                    p,
                    normal: (p - self.center) / self.radius,
                    material: self.material,
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

fn random_in_unit_sphere() -> Vector3 {
    let mut rng = rand::thread_rng();
    let mut p = Vector3::new(rng.gen::<f64>(), rng.gen::<f64>(), rng.gen::<f64>()) * 2.0
        - Vector3::new(1, 1, 1);
    while p.squared_length() >= 1.0 {
        p = Vector3::new(rng.gen::<f64>(), rng.gen::<f64>(), rng.gen::<f64>()) * 2.0
            - Vector3::new(1, 1, 1);
    }

    p
}

fn color(ray: &Ray, world: &World, depth: i64) -> Vector3 {
    if let Some(hit) = world.hit(ray, 0.001, std::f64::MAX) {
        if depth < 50 {
            if let Some((scattered, attenuation)) = hit.material.scatter(ray, &hit) {
                attenuation * color(&scattered, world, depth + 1)
            } else {
                Vector3::new(0, 0, 0)
            }
        } else {
            Vector3::new(0, 0, 0)
        }
    } else {
        let normalized_direction = ray.direction.normalized();
        let t = 0.5 * (normalized_direction.y + 1.0);

        Vector3::new(1.0, 1.0, 1.0) * (1.0 - t) + Vector3::new(0.5, 0.7, 1.0) * t
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

    let sphere1 = Sphere::new(
        Vector3::new(0.0, 0.0, -1.0),
        0.5,
        Material::Lambertian {
            albedo: Vector3::new(0.8, 0.3, 0.3),
        },
    );
    let sphere2 = Sphere::new(
        Vector3::new(0.0, -100.5, -1.0),
        100.0,
        Material::Lambertian {
            albedo: Vector3::new(0.8, 0.8, 0.0),
        },
    );
    let sphere3 = Sphere::new(
        Vector3::new(1.0, 0.0, -1.0),
        0.5,
        Material::Metal {
            albedo: Vector3::new(0.8, 0.6, 0.2),
            fuzz: 0.3,
        },
    );
    let sphere4 = Sphere::new(
        Vector3::new(-1.0, 0.0, -1.0),
        0.5,
        Material::Metal {
            albedo: Vector3::new(0.8, 0.8, 0.8),
            fuzz: 1.0,
        },
    );

    let world = World(vec![
        Box::new(sphere1),
        Box::new(sphere2),
        Box::new(sphere3),
        Box::new(sphere4),
    ]);
    let camera = Camera::new();

    for y in (0..height).rev() {
        for x in 0..width {
            let pixel = (0..num_samples).fold(Color::new(0, 0, 0), |result, _| {
                let u = (f64::from(x) + rng.gen::<f64>()) / f64::from(width);
                let v = (f64::from(y) + rng.gen::<f64>()) / f64::from(height);

                let r = camera.ray(u, v);
                let color = Color::from(color(&r, &world, 0));

                result + color
            }) / f64::from(num_samples);

            let pixel = Color::new(pixel.r.sqrt(), pixel.g.sqrt(), pixel.b.sqrt());

            println!(
                "{} {} {}",
                (pixel.r * 255.99) as u8,
                (pixel.g * 255.99) as u8,
                (pixel.b * 255.99) as u8
            );
        }
    }
}
