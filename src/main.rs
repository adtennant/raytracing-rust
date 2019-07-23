#![allow(unknown_lints)]
#![warn(clippy::all)]

use rand::prelude::*;
use rayon::prelude::*;
use std::io::{self, Write};

#[derive(Clone, Copy)]
struct Vector3 {
    x: f64,
    y: f64,
    z: f64,
}

impl Vector3 {
    const ZERO: Vector3 = Vector3 {
        x: 0.0,
        y: 0.0,
        z: 0.0,
    };
    const ONE: Vector3 = Vector3 {
        x: 1.0,
        y: 1.0,
        z: 1.0,
    };

    fn new(x: impl Into<f64>, y: impl Into<f64>, z: impl Into<f64>) -> Self {
        Vector3 {
            x: x.into(),
            y: y.into(),
            z: z.into(),
        }
    }

    fn random_in_unit_sphere() -> Self {
        let mut rng = rand::thread_rng();

        loop {
            let p = Vector3::new(rng.gen::<f64>(), rng.gen::<f64>(), rng.gen::<f64>()) * 2.0
                - Vector3::ONE;

            if p.squared_length() < 1.0 {
                break p;
            }
        }
    }

    fn cross(&self, other: &Vector3) -> Vector3 {
        Vector3::new(
            self.y * other.z - self.z * other.y,
            -(self.x * other.z - self.z * other.x),
            self.x * other.y - self.y * other.x,
        )
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

    fn refract(&self, n: &Vector3, ni_over_nt: f64) -> Option<Vector3> {
        let uv = self.normalized();
        let dt = uv.dot(n);
        let discriminant = 1.0 - ni_over_nt * ni_over_nt * (1.0 - dt * dt);

        if discriminant > 0.0 {
            Some((uv - *n * dt) * ni_over_nt - *n * discriminant.sqrt())
        } else {
            None
        }
    }

    fn squared_length(&self) -> f64 {
        self.x * self.x + self.y * self.y + self.z * self.z
    }
}

impl std::ops::Neg for Vector3 {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Vector3 {
            x: -self.x,
            y: -self.y,
            z: -self.z,
        }
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

fn schlick(cosine: f64, refractive_index: f64) -> f64 {
    let r0 = (1.0 - refractive_index) / (1.0 + refractive_index);
    let r0 = r0 * r0;

    r0 + (1.0 - r0) * (1.0 - cosine).powf(5.0)
}

#[derive(Clone, Copy)]
enum Material {
    Lambertian { albedo: Vector3 },
    Metal { albedo: Vector3, fuzz: f64 },
    Dielectric { refractive_index: f64 },
}

impl Material {
    fn lambertian(albedo: Vector3) -> Self {
        Material::Lambertian { albedo }
    }

    fn metal(albedo: Vector3, fuzz: impl Into<f64>) -> Self {
        Material::Metal {
            albedo,
            fuzz: fuzz.into(),
        }
    }

    fn dielectric(refractive_index: impl Into<f64>) -> Self {
        Material::Dielectric {
            refractive_index: refractive_index.into(),
        }
    }

    fn scatter(&self, ray: &Ray, hit: &Hit) -> Option<(Ray, Vector3)> {
        match *self {
            Material::Lambertian { albedo } => {
                let target = hit.position + hit.normal + Vector3::random_in_unit_sphere();
                let scattered = Ray::new(hit.position, target - hit.position);
                Some((scattered, albedo))
            }
            Material::Metal { albedo, fuzz } => {
                let reflected = ray.direction.normalized().reflect(&hit.normal);
                let scattered = Ray::new(
                    hit.position,
                    reflected + Vector3::random_in_unit_sphere() * fuzz,
                );

                if scattered.direction.dot(&hit.normal) > 0.0 {
                    Some((scattered, albedo))
                } else {
                    None
                }
            }
            Material::Dielectric { refractive_index } => {
                let (outward_normal, ni_over_nt, cosine) = if ray.direction.dot(&hit.normal) > 0.0 {
                    (
                        -hit.normal,
                        refractive_index,
                        refractive_index * ray.direction.dot(&hit.normal) / ray.direction.length(),
                    )
                } else {
                    (
                        hit.normal,
                        1.0 / refractive_index,
                        -ray.direction.dot(&hit.normal) / ray.direction.length(),
                    )
                };

                let scattered = match ray.direction.refract(&outward_normal, ni_over_nt) {
                    Some(refracted) if random::<f64>() >= schlick(cosine, refractive_index) => {
                        Ray::new(hit.position, refracted)
                    }
                    _ => Ray::new(hit.position, ray.direction.reflect(&hit.normal)),
                };
                let attenuation = Vector3::ONE;

                Some((scattered, attenuation))
            }
        }
    }
}

struct Hit {
    distance: f64,
    position: Vector3,
    normal: Vector3,
    material: Material,
}

enum Shape {
    Sphere {
        center: Vector3,
        radius: f64,
        material: Material,
    },
}

impl Shape {
    fn sphere(center: Vector3, radius: impl Into<f64>, material: Material) -> Self {
        Shape::Sphere {
            center,
            radius: radius.into(),
            material,
        }
    }

    fn hit(
        &self,
        ray: &Ray,
        min_distance: impl Into<f64>,
        max_distance: impl Into<f64>,
    ) -> Option<Hit> {
        let min_distance = min_distance.into();
        let max_distance = max_distance.into();

        match *self {
            Shape::Sphere {
                center,
                radius,
                material,
            } => {
                let oc = ray.origin - center;

                let a = ray.direction.dot(&ray.direction);
                let b = oc.dot(&ray.direction);
                let c = oc.dot(&oc) - radius * radius;

                let discriminant = b * b - a * c;

                if discriminant > 0.0 {
                    let distance = (-b - discriminant.sqrt()) / a;

                    if distance > min_distance && distance < max_distance {
                        let position = ray.point_at(distance);

                        return Some(Hit {
                            distance,
                            position,
                            normal: (position - center) / radius,
                            material,
                        });
                    }

                    let distance = (-b + discriminant.sqrt()) / a;

                    if distance > min_distance && distance < max_distance {
                        let position = ray.point_at(distance);

                        return Some(Hit {
                            distance,
                            position,
                            normal: (position - center) / radius,
                            material,
                        });
                    }
                }

                None
            }
        }
    }
}

struct World(Vec<Shape>);

impl World {
    fn hit(
        &self,
        ray: &Ray,
        min_distance: impl Into<f64>,
        max_distance: impl Into<f64>,
    ) -> Option<Hit> {
        let min_distance = min_distance.into();
        let max_distance = max_distance.into();

        let mut nearest = max_distance;
        self.0
            .iter()
            .filter_map(|h| h.hit(ray, min_distance, max_distance))
            .fold(None, |result, h| {
                if h.distance < nearest {
                    nearest = h.distance;
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
    right: Vector3,
    up: Vector3,
    lens_radius: f64,
}

impl Camera {
    fn new(
        origin: Vector3,
        target: Vector3,
        up: Vector3,
        vfov: impl Into<f64>,
        aspect: impl Into<f64>,
        aperture: impl Into<f64>,
        focus_distance: impl Into<f64>,
    ) -> Self {
        let theta = vfov.into() * std::f64::consts::PI / 180.0;
        let half_height = (theta / 2.0).tan();
        let half_width = aspect.into() * half_height;

        let forward = (origin - target).normalized();
        let right = up.cross(&forward).normalized();
        let up = forward.cross(&right);

        let focus_distance = focus_distance.into();

        Camera {
            origin,
            lower_left: origin
                - right * half_width * focus_distance
                - up * half_height * focus_distance
                - forward * focus_distance,
            horizontal: right * 2.0 * half_width * focus_distance,
            vertical: up * 2.0 * half_height * focus_distance,
            right,
            up,
            lens_radius: aperture.into() / 2.0,
        }
    }

    fn ray(&self, u: impl Into<f64>, v: impl Into<f64>) -> Ray {
        let rd = Vector3::random_in_unit_sphere() * self.lens_radius;
        let offset = self.right * rd.x + self.up * rd.y;

        Ray::new(
            self.origin + offset,
            self.lower_left + self.horizontal * u.into() + self.vertical * v.into()
                - self.origin
                - offset,
        )
    }
}

fn color(ray: &Ray, world: &World, depth: i64) -> Vector3 {
    if let Some(hit) = world.hit(ray, 0.001, std::f64::MAX) {
        match hit.material.scatter(ray, &hit) {
            Some((ref scattered, attenuation)) if depth < 50 => {
                attenuation * color(scattered, world, depth + 1)
            }
            _ => Vector3::ZERO,
        }
    } else {
        let normalized_direction = ray.direction.normalized();
        let t = 0.5 * (normalized_direction.y + 1.0);

        Vector3::ONE * (1.0 - t) + Vector3::new(0.5, 0.7, 1) * t
    }
}

fn random_scene() -> World {
    let mut scene = vec![
        Shape::sphere(
            Vector3::new(0, -1000, -1),
            1000,
            Material::lambertian(Vector3::new(0.5, 0.5, 0.5)),
        ),
        Shape::sphere(Vector3::new(0, 1, 0), 1, Material::dielectric(1.5)),
        Shape::sphere(
            Vector3::new(-4, 1, 0),
            1,
            Material::lambertian(Vector3::new(0.4, 0.2, 0.1)),
        ),
        Shape::sphere(
            Vector3::new(4, 1, 0),
            1,
            Material::metal(Vector3::new(0.7, 0.6, 0.5), 0),
        ),
    ];

    let mut rng = rand::thread_rng();

    for a in -11..11 {
        for b in -11..11 {
            let center = Vector3::new(
                f64::from(a) + 0.9 * rng.gen::<f64>(),
                0.2,
                f64::from(b) + 0.9 + rng.gen::<f64>(),
            );

            if (center - Vector3::new(4, 0.2, 0)).length() > 0.9 {
                let material = match rng.gen::<f64>() {
                    x if x >= 0.0 && x < 0.8 => Material::lambertian(Vector3::new(
                        rng.gen::<f64>() * rng.gen::<f64>(),
                        rng.gen::<f64>() * rng.gen::<f64>(),
                        rng.gen::<f64>() * rng.gen::<f64>(),
                    )),
                    x if x >= 0.8 && x < 0.95 => Material::metal(
                        Vector3::new(
                            0.5 * (1.0 + rng.gen::<f64>()),
                            0.5 * (1.0 + rng.gen::<f64>()),
                            0.5 * (1.0 + rng.gen::<f64>()),
                        ),
                        0.5 * rng.gen::<f64>(),
                    ),
                    _ => Material::dielectric(1.5),
                };

                scene.push(Shape::sphere(center, 0.2, material));
            }
        }
    }

    World(scene)
}

fn main() -> io::Result<()> {
    let width = 200;
    let height = 100;
    let num_samples = 100;

    let world = random_scene();

    let start = std::time::Instant::now();

    let origin = Vector3::new(16, 2, 4);
    let target = Vector3::ZERO;
    let camera = Camera::new(
        origin,
        target,
        Vector3::new(0, 1, 0),
        15,
        f64::from(width) / f64::from(height),
        0.2,
        (origin - target).length(),
    );

    let header = format!("P3\n{} {}\n255", width, height);
    let pixels = (0..height * width)
        .into_par_iter()
        .map(|i| {
            let mut rng = rand::thread_rng();

            let y = height - i / width - 1;
            let x = i % width;

            // no noticable performance different switching to par_iter here, likely due to the fold
            let pixel = (0..num_samples).fold(Color::new(0, 0, 0), |result, _| {
                let u = (f64::from(x) + rng.gen::<f64>()) / f64::from(width);
                let v = (f64::from(y) + rng.gen::<f64>()) / f64::from(height);

                let ray = camera.ray(u, v);
                let color = Color::from(color(&ray, &world, 0));

                result + color
            }) / f64::from(num_samples);

            let pixel = Color::new(pixel.r.sqrt(), pixel.g.sqrt(), pixel.b.sqrt());

            format!(
                "{} {} {}",
                (pixel.r * 255.0) as u8,
                (pixel.g * 255.0) as u8,
                (pixel.b * 255.0) as u8
            )
        })
        .collect::<Vec<_>>()
        .join("\n");

    let ppm = format!("{}\n{}\n", header, pixels);
    io::stdout().write_all(ppm.as_bytes())?;

    let end = std::time::Instant::now();
    io::stderr().write_all(format!("Finished rendering in {:?}.\n", end - start).as_bytes())?;

    Ok(())
}
