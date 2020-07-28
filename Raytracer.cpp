// This code is cut from the rest of the file to avoid uploading starter code
// Raytracer- Anisha Braganza

// Render
hdr_image scene::render() const noexcept {

  assert(camera_);
  assert(viewport_);
  assert(projection_);
  assert(shader_);

  assert(viewport_->x_resolution() > 0);
  assert(viewport_->y_resolution() > 0);

  size_t w = viewport_->x_resolution(),
         h = viewport_->y_resolution();


  hdr_image result(w, h, background_);
  assert(!result.is_empty());

  for (size_t y = 0; y < h; ++y) {
    for (size_t x = 0; x < w; ++x) {

      // - Compute the (u, v) corresponding to (x, y)
      vector2<double> uv_values = viewport_->uv(x, y);
      // - Use the projection object to create the view ray based on that (u, v)
      view_ray current_view_ray = projection_->compute_view_ray(camera(), uv_values[0], uv_values[1]);
      // - Use the scene object to trace the view ray and find an
      //   intersection.
      // - If there is no intersection, paint result.pixel(x, y) with
      //   the scene's background color
      std::optional<intersection> front_object = intersect(current_view_ray);
      if (front_object == std::nullopt){
        result.pixel(x,y, background());
      }

      // - Otherwise, use the shader object to compute the color for
      //   result.pixel(x, y) based on the intersection object.
      else{
        result.pixel(x, y, shader_->shade(*this, camera(), *front_object));
      }
    }
  }
  return result;
}


// calculate intersection
std::optional<intersection> scene::intersect(const view_ray& ray) const noexcept {

  double t_min = 0.0;
  double t_upper_bound = std::numeric_limits<double>::infinity();

  std::optional<intersection> closest = std::nullopt;

  //exhaustive for closest object
  for (auto itr = std::begin(objects()); itr != std::end(objects()); itr++) {
    std::optional<intersection> object_hit = (*itr)->intersect(ray, t_min, t_upper_bound);
    if (object_hit.has_value()){
      double t = object_hit->t();
      if ( t >= t_min && t < t_upper_bound){
        closest = object_hit;
        t_upper_bound = t;
      }
    }
  }

  return closest;
}

constexpr camera::camera(const vector3<double>& eye,
	                       const vector3<double>& view_direction,
	                       const vector3<double>& up) noexcept {

  //compute basis, then normalize

  eye_ = eye;
  w_ = (view_direction * -1).normalized();
  u_ = (up.cross(w_)).normalized();
  v_ = (w_.cross(u_)).normalized();
}

vector2<double> viewport::uv(size_t x, size_t y) const noexcept {


  vector2<double> orthPt;

  orthPt[0] = left_ + ((right_ - left_) * (x + 0.5)) / x_resolution_;
  orthPt[1] = bottom_ + ((top_ - bottom_) * (y + 0.5)) / y_resolution_;

  return orthPt;
}

view_ray orthographic_projection::compute_view_ray(const camera& c,
					                                         double u,
					                                         double v) const noexcept {


  gfx::vector3<double> orthDirection = c.w() * -1;  // ADD NEG
  gfx::vector3<double> orthOrigin = c.eye() + (c.u()*u) + (c.v()*v);

  return view_ray(orthOrigin, orthDirection);
}

view_ray perspective_projection::compute_view_ray(const camera& c,
					                                        double u,
					                                        double v) const noexcept {


  gfx::vector3<double> perOrigin = c.eye();
  // -dw + uU + vV
  gfx::vector3<double> perDirection = c.w()*(focal_length() * -1) + c.u()*u + c.v()*v; //ADD NEG


  return view_ray(perOrigin, perDirection);
}

hdr_rgb flat_shader::shade(const scene& scene,
			                     const camera& camera,
			                     const intersection& xsect) const noexcept {

  return xsect.object().color();
}


hdr_rgb blinn_phong_shader::shade(const scene& scene,
				                          const camera& camera,
				                          const intersection& xsect) const noexcept {


  gfx::hdr_rgb ambientLight;
  double diffuseLight = 0.0, specularLight = 0.0;

  ambientLight.r(ambient_coefficient()* ambient_color().r());
  ambientLight.b(ambient_coefficient()* ambient_color().b());
  ambientLight.g(ambient_coefficient()* ambient_color().g());

  // interate through scene lights
  for (auto itr = std::begin(scene.lights()); itr != std::end(scene.lights()); itr++){
    // construct l(i) then normalize
    vector3<double> lightLocationVector = ((*itr)->location() - xsect.location()).normalized();
    // construct h(i)
    vector3<double> bisector = (lightLocationVector + (camera.eye() - xsect.location()).normalized()).normalized();

    // diffuse light
    diffuseLight += diffuse_coefficient() * std::max(0.0, xsect.normal() * lightLocationVector);
    //specular Light
    specularLight += specular_coefficient() * std::pow(std::max(0.0, xsect.normal() * bisector), xsect.object().shininess());
  }

  // calculate rgb values
  float r =  (ambientLight.r() + specularLight + xsect.object().color().r() *diffuseLight);// + xsect.object().color().r() ;
  float g =(ambientLight.g() + specularLight +  xsect.object().color().g() *diffuseLight);// + xsect.object().color().g() ;
  float b = (ambientLight.b() + specularLight + xsect.object().color().b() *diffuseLight);// + xsect.object().color().b() ;

  // bounds tests [0, 1]
  r = (r > 1.0f) ? 1.0f : std::max(0.0f, r),
  g = (g > 1.0f) ? 1.0f : std::max(0.0f, g),
  b = (b > 1.0f) ? 1.0f : std::max(0.0f, b);

  return gfx::hdr_rgb(r , g, b);
}

std::optional<intersection>
  scene_sphere::intersect(
    const view_ray& ray,
    double t_min,
    double t_upper_bound) const noexcept {

  assert(t_min < t_upper_bound);

  // e is origin()
  // d is direction()
  gfx::vector3<double> d = ray.direction();
  gfx::vector3<double> e = ray.origin();
  gfx::vector3<double> eMinC = e - center();
  double R = radius();

  // B^2 - 4AC
  double Bsqur = std::pow((d * 2) * (eMinC), 2);
  double frAC = (d * d) * (eMinC * eMinC - (R * R)) *  4;
  double discriminant = Bsqur - frAC;

  // if negative, imaginary so no intersection
  if (discriminant < 0.0) {return std::nullopt; }

  // if 0, one intersection
  else if (discriminant  == 0.0){
    // t = -d dot (e - c) / 2a
    double t = ((d * -1.0) * (e - center()))  / (d * d);
    vector3<double> p = ray.origin() + d * t;
    vector3<double> n = (p - center()) * 2.0;
    return gfx::intersection(this, p, n.normalized(), t);
  }
  // if postiive, 2 solutions
  else if (discriminant > 0.0){
    double t1 = (((d * -1)* (e - center())) + std::sqrt(discriminant)) / (d * d);
    double t2 = (((d * -1)* (e - center())) - std::sqrt(discriminant)) / (d * d);

    double t = (t1 < t2) ? t1 : t2;  // find first intersectoin

    vector3<double> p = ray.origin() + d * t;
    vector3<double> n = (p - center()) * 2;

    return gfx::intersection(this, p, n.normalized(), t);
  }

  return std::nullopt;
}

std::optional<intersection>
  scene_triangle::intersect(
    const view_ray& ray,
    double t_min,
    double t_upper_bound) const noexcept {

  assert(t_min < t_upper_bound);

  // vars
  vector3<double> a = this->a();
  vector3<double> b = this->b();
  vector3<double> c = this->c();
  vector3<double> d = ray.direction();
  vector3<double> e = ray.origin();

  gfx::matrix<double, 3, 3> matrixA {a[0] - b[0], a[0] - c[0], d[0],
                                     a[1] - b[1], a[1] - c[1], d[1],
                                     a[2] - b[2], a[2] - c[2], d[2]};

  vector3<double> solution = {a[0] - e[0], a[1] - e[1], a[2] - e[2]};
  vector3<double> solutions = matrixA.solve(solution);

  // bound intersection checking
  double t = solutions[2];
  if ((t < t_min) || (t > t_upper_bound)) { return std::nullopt;}

  double gamma = solutions[1];
  if ((gamma < 0) || (gamma > 1)) { return std::nullopt;}

  double beta = solutions[0];
  if ((beta < 0) || (beta > (1 - gamma))) { return std::nullopt;}

  vector3<double> p = (ray.origin() + d * t);
  vector3<double> n = ((b-a).cross(c-a));

  return  gfx::intersection(this, p, n.normalized(), t);;
}
