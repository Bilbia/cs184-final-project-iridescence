#include "bsdf.h"

#include <algorithm>
#include <iostream>
#include <utility>

#include "application/visual_debugger.h"

using std::max;
using std::min;
using std::swap;

namespace CGL {

// Mirror BSDF //

    Vector3D MirrorBSDF::f(const Vector3D wo, const Vector3D wi) {
        return Vector3D();
    }

    Vector3D MirrorBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {

        // TODO Project 3-2: Part 1
        // Implement MirrorBSDF
//        *pdf = 1;
//        reflect(wo, wi);
//        return reflectance/abs_cos_theta(*wi);
      if (!refract(wo, wi, 1)) {
        reflect(wo, wi);
        *pdf = 1;
        return reflectance / abs_cos_theta(*wi);
      }
      double R0 = pow((1 - 1) / (1 + 1), 2);
      double R = R0 + (1 - R0) * pow(1 - abs_cos_theta(wo), 5);
      if (coin_flip(R)) {
        reflect(wo, wi);
        *pdf = R;
        return R * reflectance / abs_cos_theta(*wi);
      }
      double eta = 1;
      if (wo.z >= 0) {
        eta = 1. / 1;
      }
      refract(wo, wi, 1);
      *pdf = 1 - R;
      Vector3D transmittance = (1,1,1);
      return (1 - R) * transmittance / abs_cos_theta(*wi) / pow(eta, 2);



    }

    void MirrorBSDF::render_debugger_node()
    {
        if (ImGui::TreeNode(this, "Mirror BSDF"))
        {
            DragDouble3("Reflectance", &reflectance[0], 0.005);
            ImGui::TreePop();
        }
    }

// Microfacet BSDF //

    double MicrofacetBSDF::G(const Vector3D wo, const Vector3D wi) {
        return 1.0 / (1.0 + Lambda(wi) + Lambda(wo));
    }

    double MicrofacetBSDF::D(const Vector3D h) {
        // TODO Project 3-2: Part 2
        // Compute Beckmann normal distribution function (NDF) here.
        // You will need the roughness alpha.

        //double theta = getTheta(h);
        //return exp(-(pow(tan(theta), 2) / pow(alpha, 2))) / (PI * pow(alpha, 2) * pow(cos(theta), 4));
        double theta_h = acos(h.z);
        return  exp(-1. * pow(tan(theta_h) / alpha, 2)) / (PI * pow(alpha, 2) * pow(h.z, 4));
    }

    // Helper functions for Iridescence. R and T functions represent reflection and transmission coefficients
    Vector3D costrans(Vector3D sini) {
        double mag = std::sqrt(sini.x * sini.x + sini.y * sini.y + sini.z * sini.z);
        Vector3D cos;
        cos.x = std::sqrt(1 - (sini.x / mag) * (sini.x / mag));
        cos.y = std::sqrt(1 - (sini.y / mag) * (sini.y / mag));
        cos.z = std::sqrt(1 - (sini.z / mag) * (sini.z / mag));
        return cos;
    }

    Vector3D sintrans(Vector3D eta, Vector3D eta2, Vector3D cos) {
        Vector3D sin = (eta / eta2) * (eta / eta2);

        sin.x = sin.x * (1 - cos.x * cos.x);
        sin.y = sin.y * (1 - cos.y * cos.y);
        sin.z = sin.z * (1 - cos.z * cos.z);
        return sin;
    }

    Vector3D Rs(Vector3D eta, Vector3D eta2, Vector3D cosi, Vector3D cost) {
        return (eta * cosi - eta2 * cost) / (eta * cosi + eta2 * cost);
    }

    Vector3D Rp(Vector3D eta, Vector3D eta2, Vector3D cosi, Vector3D cost) {
        return (eta2 * cosi - eta * cost) / (eta2 * cosi + eta * cost);
    }

    Vector3D Ts(Vector3D eta, Vector3D eta2, Vector3D cosi, Vector3D cost) {
        return (2 * eta * cosi) / (eta * cosi + eta2 * cost);
    }

    Vector3D Tp(Vector3D eta, Vector3D eta2, Vector3D cosi, Vector3D cost) {
        return (2 * eta * cosi) / (eta2 * cosi + eta * cost);
    }

    Vector3D Psi(Vector3D lambda, double thick, Vector3D eta0, Vector3D eta1, Vector3D eta2, Vector3D costheta) {
        // 180 phase shift if eta2>eta
        Vector3D delta1 = 0;
        if (eta1.x > eta0.x) delta1.x = PI;
        if (eta1.y > eta0.y) delta1.y = PI;
        if (eta1.z > eta0.z) delta1.z = PI;

        Vector3D delta2 = 0;
        if (eta1.x > eta2.x) delta2.x = PI;
        if (eta1.y > eta2.y) delta2.y = PI;
        if (eta1.z > eta2.z) delta2.z = PI;

        Vector3D delta = delta1 + delta2;
        return (2 * PI / lambda) * (2 * eta1 * thick * costheta) + delta;
    }


    Vector3D MicrofacetBSDF::F(const Vector3D wi) {
        // TODO Project 3-2: Part 2
        // Compute Fresnel term for reflection on dielectric-conductor interface.
        // You will need both eta and etaK, both of which are Vector3D.

        // IRIDESCENCE
        //std::cout << eta << std::endl;

        // Fresnell R and T
        Vector3D eta0 = 1;
        Vector3D eta2 = 1.33;
        Vector3D cos0 = wi.z;
        Vector3D sin1 = sintrans(eta0, eta, cos0);
        Vector3D cos1 = costrans(sin1);
        Vector3D sin2 = sintrans(eta, eta2, cos1);
        Vector3D cos2 = costrans(sin2);

        Vector3D alphaS = Rs(eta, eta0, cos1, cos0) * Rs(eta, eta2, cos1, cos2);
        Vector3D alphaP = Rp(eta, eta0, cos1, cos0) * Rp(eta, eta2, cos1, cos2);
        Vector3D betaS = Ts(eta0, eta, cos0, cos1) * Ts(eta, eta2, cos1, cos2);
        Vector3D betaP = Tp(eta0, eta, cos0, cos1) * Tp(eta, eta2, cos1, cos2);

        // Phase changes
        Vector3D lambda = { 614, 549, 466 }; //wavelengths of R,G,B

        double thick = 200 + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / (400.0 - 300.0))); //random thickness between 200 and 400
        Vector3D psi = Psi(lambda, thick, eta0, eta, eta2, cos1);
        Vector3D cosPsi = { cos(psi.x), cos(psi.y), cos(psi.z) };
        //std::cout << cosPsi << std::endl;

        Vector3D ts = ((betaS * betaS) / (alphaS * alphaS - 2 * alphaS * cosPsi + 1));
        Vector3D tp = ((betaP * betaP) / (alphaP * alphaP - 2 * alphaP * cosPsi + 1));
        Vector3D T = ((eta2 * cos2) / (eta0 * cos0)) * (ts + tp) / 2;
        return Vector3D(1.0) - T;
    }

    Vector3D MicrofacetBSDF::f(const Vector3D wo, const Vector3D wi) {
        // TODO Project 3-2: Part 2
        // Implement microfacet model here.
        //if (wi.z > 0 && wo.z > 0) {
        //    Vector3D h = (wo + wi) / (wo + wi).norm();
        //    return (F(wi) * G(wo, wi) * D(h)) / (4 * wo.z * wi.z);
        //}

        //return Vector3D();
        if (wo.z > 0 && wi.z > 0) {
            Vector3D h = (wo + wi) / (wo + wi).norm();
            Vector3D f = (F(wi) * G(wo, wi) * D(h)) / (4. * wo.z * wi.z);
            return f;
        }
        return Vector3D(0.);
    }

    Vector3D MicrofacetBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {
        // TODO Project 3-2: Part 2
        // *Importance* sample Beckmann normal distribution function (NDF) here.
        // Note: You should fill in the sampled direction *wi and the corresponding *pdf,
        //       and return the sampled BRDF value.

        Vector2D r = sampler.get_sample();
        double theta = atan(sqrt(-pow(alpha, 2) * log(1. - r.x)));
        double phi = 2. * PI * r.y;
        Vector3D h = Vector3D(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta));
        double ptheta = (2. * sin(theta) / (pow(alpha, 2) * pow(cos(alpha), 3))) * exp(-pow(tan(theta), 2) / pow(alpha, 2));
        double pphi = 1. / (2. * PI);
        *wi = -wo + 2. * dot(wo, h) * h;
        wi->normalize();
        *pdf = ((ptheta * pphi) / sin(theta))/(4. * dot(*wi, h));
        return  MicrofacetBSDF::f(wo, *wi);
    }

    void MicrofacetBSDF::render_debugger_node()
    {
        if (ImGui::TreeNode(this, "Micofacet BSDF"))
        {
            DragDouble3("eta", &eta[0], 0.005);
            DragDouble3("K", &k[0], 0.005);
            DragDouble("alpha", &alpha, 0.005);
            ImGui::TreePop();
        }
    }

// Refraction BSDF //

    Vector3D RefractionBSDF::f(const Vector3D wo, const Vector3D wi) {
        return Vector3D();
    }

    Vector3D RefractionBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {
        // TODO Project 3-2: Part 1
        // Implement RefractionBSDF

        if (refract(wo, wi, ior)) {
            *pdf = 1;
            double eta = ior;
            if (wo.z >= 0) {
                eta = 1. / ior;
            }
            return transmittance / abs_cos_theta(*wi) / pow(eta, 2);
        }
        return Vector3D();
    }

    void RefractionBSDF::render_debugger_node()
    {
        if (ImGui::TreeNode(this, "Refraction BSDF"))
        {
            DragDouble3("Transmittance", &transmittance[0], 0.005);
            DragDouble("ior", &ior, 0.005);
            ImGui::TreePop();
        }
    }

// Glass BSDF //

    Vector3D GlassBSDF::f(const Vector3D wo, const Vector3D wi) {
        return Vector3D();
    }

    Vector3D GlassBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {

        // TODO Project 3-2: Part 1
        // Compute Fresnel coefficient and either reflect or refract based on it.

        // compute Fresnel coefficient and use it as the probability of reflection
        // - Fundamentals of Computer Graphics page 305
        if (!refract(wo, wi, 1)) {
            reflect(wo, wi);
            *pdf = 1;
            return reflectance / abs_cos_theta(*wi);
        }
        double R0 = pow((1 - 1) / (1 + 1), 2);
        double R = R0 + (1 - R0) * pow(1 - abs_cos_theta(wo), 5);
        if (coin_flip(R)) {
            reflect(wo, wi);
            *pdf = R;
            return R * reflectance / abs_cos_theta(*wi);
        }
        double eta = 1;
        if (wo.z >= 0) {
            eta = 1. / 1;
        }
        refract(wo, wi, 1);
        *pdf = 1 - R;
        Vector3D transmittance = (1,1,1);
        return (1 - R) * transmittance / abs_cos_theta(*wi) / pow(eta, 2);
//        *pdf = 1;
//        reflect(wo, wi);
//        return reflectance/abs_cos_theta(*wi);

    }

    void GlassBSDF::render_debugger_node()
    {
        if (ImGui::TreeNode(this, "Refraction BSDF"))
        {
            DragDouble3("Reflectance", &reflectance[0], 0.005);
            DragDouble3("Transmittance", &transmittance[0], 0.005);
            DragDouble("ior", &ior, 0.005);
            ImGui::TreePop();
        }
    }

    void BSDF::reflect(const Vector3D wo, Vector3D* wi) {

        // TODO Project 3-2: Part 1
        // Implement reflection of wo about normal (0,0,1) and store result in wi.
        *wi = wo * Vector3D(-1, -1, 1);

    }

    bool BSDF::refract(const Vector3D wo, Vector3D* wi, double ior) {

        // TODO Project 3-2: Part 1
        // Use Snell's Law to refract wo surface and store result ray in wi.
        // Return false if refraction does not occur due to total internal reflection
        // and true otherwise. When dot(wo,n) is positive, then wo corresponds to a
        // ray entering the surface through vacuum.

        double eta = ior;
        if (wo.z >= 0) {
            eta = 1. / ior;
        }

        double sn = 1 - (pow(eta, 2) * (1 - pow(wo.z, 2)));
        if (sn < 0) {
            return false;
        }

        wi->x = -eta * wo.x;
        wi->y = -eta * wo.y;
        wi->z = sqrt(sn);
        if (wo.z >= 0) wi->z *= -1;

        return true;
    }

} // namespace CGL
