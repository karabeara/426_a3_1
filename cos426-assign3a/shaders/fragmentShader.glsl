// set the precision of the float values (necessary if using float)
#ifdef GL_FRAGMENT_PRECISION_HIGH
precision highp float;
#else
precision mediump float;
#endif
precision mediump int;

// define constant parameters
// EPS is for the precision issue (see precept slide)
#define INFINITY 1.0e+12
#define EPS 1.0e-3

// define constants for scene setting
#define MAX_LIGHTS 10

// define texture types
#define NONE 0
#define CHECKERBOARD 1
#define MYSPECIAL 2

// define material types
#define BASICMATERIAL 1
#define PHONGMATERIAL 2
#define LAMBERTMATERIAL 3

// define reflect types - how to bounce rays
#define NONEREFLECT 1
#define MIRRORREFLECT 2
#define GLASSREFLECT 3

struct Shape {
    int shapeType;
    vec3 v1;
    vec3 v2;
    float rad;
};

struct Material {
    int materialType;
    vec3 color;
    float shininess;
    vec3 specular;

    int materialReflectType;
    float reflectivity;
    float refractionRatio;
    int special;

};

struct Object {
    Shape shape;
    Material material;
};

struct Light {
    vec3 position;
    vec3 color;
    float intensity;
    float attenuate;
};

struct Ray {
    vec3 origin;
    vec3 direction;
};

struct Intersection {
    vec3 position;
    vec3 normal;
};

// uniform
uniform mat4 uMVMatrix;
uniform int frame;
uniform float height;
uniform float width;
uniform vec3 camera;
uniform int numObjects;
uniform int numLights;
uniform Light lights[MAX_LIGHTS];
uniform vec3 objectNorm;

varying vec2 v_position;

// find then position some distance along a ray
vec3 rayGetOffset( Ray ray, float dist ) {
    return ray.origin + ( dist * ray.direction );
}

// if a newly found intersection is closer than the best found so far, record the new intersection and return true;
// otherwise leave the best as it was and return false.
bool chooseCloserIntersection( float dist, inout float best_dist, inout Intersection intersect, inout Intersection best_intersect ) {
    if ( best_dist <= dist ) return false;
    best_dist = dist;
    best_intersect.position = intersect.position;
    best_intersect.normal   = intersect.normal;
    return true;
}

// put any general convenience functions you want up here
// ----------- STUDENT CODE BEGIN ------------
// ----------- Our reference solution uses 135 lines of code.
// ----------- STUDENT CODE END ------------

// forward declaration
float rayIntersectScene( Ray ray, out Material out_mat, out Intersection out_intersect );

// Plane
// this function can be used for plane, triangle, and box
float findIntersectionWithPlane( Ray ray, vec3 norm, float dist, out Intersection intersect ) {
    float a   = dot( ray.direction, norm );
    float b   = dot( ray.origin, norm ) - dist;

    if ( a < 0.0 && a > 0.0 ) return INFINITY;

    float len = -b/a;
    if ( len < EPS ) return INFINITY;

    intersect.position = rayGetOffset( ray, len );
    intersect.normal   = norm;
    return len;
}

float areaOfTriangle(vec3 v1, vec3 v2, vec3 v3, vec3 norm) {
  vec3 cross1 = cross(v2 - v1, v3 - v1);
  float area1 = 0.5 * dot(cross1, norm);
  return area1;
}

// Triangle
float findIntersectionWithTriangle( Ray ray, vec3 t1, vec3 t2, vec3 t3, out Intersection intersect ) {
    // ----------- STUDENT CODE BEGIN ------------

    // Find P - intersection of ray with point on triangle's plane
    // find normal of given triangle face
    vec3 triangleVec1 = t2 - t1;
    vec3 triangleVec2 = t3 - t1;
    vec3 triangleNormal = cross(triangleVec1, triangleVec2);
    vec3 normalized_triangleNormal = normalize(triangleNormal);

    float dist = dot(normalized_triangleNormal, t1);
    //dist = abs(dist);

    float distToPlane = findIntersectionWithPlane(ray, normalized_triangleNormal, dist, intersect);

    vec3 P = intersect.position;

    // Area of triangle
    float totalTriArea = areaOfTriangle(t1, t2, t3, normalized_triangleNormal);

    // Area of subtriangle 1
    float area1 = areaOfTriangle(t1, t2, P, normalized_triangleNormal) / totalTriArea;

    // Area of subtriangle 2
    float area2 = areaOfTriangle(t1, P, t3, normalized_triangleNormal) / totalTriArea;

    //float areaSum = abs(a) + abs(b) + abs(c);
    if (area1 >= 0.0 && area1 <= 1.0 && area2 >= 0.0 && area2 <= 1.0 && area1 + area2 >= 0.0 && area1 + area2 <= 1.0) { return distToPlane; }

    return INFINITY;

    // ----------- Our reference solution uses 22 lines of code.
    // return INFINITY; // currently reports no intersection
    // ----------- STUDENT CODE END ------------
}

// Sphere
float findIntersectionWithSphere( Ray ray, vec3 center, float radius, out Intersection intersect ) {
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 23 lines of code.
    vec3 lengthToCenter = center - ray.origin;
    vec3 normalizedDirection = normalize(ray.direction);
    float tCA = dot(lengthToCenter, normalizedDirection) ;

    if (tCA < 0.0) { return INFINITY; }
    float rSquared = radius * radius;
    float dSquared = dot(lengthToCenter,lengthToCenter) - (tCA * tCA);
    if ( dSquared > rSquared ) { return INFINITY; }

    float tHC = sqrt(rSquared - dSquared);
    float t1 = tCA - tHC;
    float t2 = tCA + tHC;
    float t;
    if (t1 < t2 && t1 > 0.0) { t = t1; }
    else                     { t = t2; }
    intersect.position = ray.origin + (t * normalizedDirection);
    intersect.normal = normalize(intersect.position - center);

    return length(t * normalizedDirection);
    // ----------- STUDENT CODE END ------------
}

// returns an array of four vec3's which are the faces 

// takes array of faces and then calculates normal and distance,
// get these two to get intersections

// check if points are on box 
// compare nomals


// Box
float findIntersectionWithBox( Ray ray, vec3 pmin, vec3 pmax, out Intersection out_intersect ) {
    // ----------- STUDENT CODE BEGIN ------------
    // pmin and pmax represent two bounding points of the box
    // pmin stores [xmin, ymin, zmin] and pmax stores [xmax, ymax, zmax]


    function()
    /* function that takes 3 points and then computes normal and then checks 
    intersection and then checks if point is on plane. */



    // get all of the faces, calculate normals , and then return the smallest one. 



    // ----------- Our reference solution uses 24 lines of code.
    return INFINITY; // currently reports no intersection
    // ----------- STUDENT CODE END ------------
}

// Cylinder
float getIntersectOpenCylinder( Ray ray, vec3 center, vec3 axis, float len, float rad, out Intersection intersect ) {
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 31 lines of code.
    return INFINITY; // currently reports no intersection
    // ----------- STUDENT CODE END ------------
}

float getIntersectDisc( Ray ray, vec3 center, vec3 norm, float rad, out Intersection intersect ) {
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 15 lines of code.
    return INFINITY; // currently reports no intersection
    // ----------- STUDENT CODE END ------------
}


float findIntersectionWithCylinder( Ray ray, vec3 center, vec3 apex, float radius, out Intersection out_intersect ) {
    vec3 axis = apex - center;
    float len = length( axis );
    axis = normalize( axis );

    Intersection intersect;
    float best_dist = INFINITY;
    float dist;

    // -- infinite cylinder
    dist = getIntersectOpenCylinder( ray, center, axis, len, radius, intersect );
    chooseCloserIntersection( dist, best_dist, intersect, out_intersect );

    // -- two caps
    dist = getIntersectDisc( ray, center, axis, radius, intersect );
    chooseCloserIntersection( dist, best_dist, intersect, out_intersect );
    dist = getIntersectDisc( ray,   apex, axis, radius, intersect );
    chooseCloserIntersection( dist, best_dist, intersect, out_intersect );

    return best_dist;
}

// Cone
float getIntersectOpenCone( Ray ray, vec3 apex, vec3 axis, float len, float radius, out Intersection intersect ) {
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 31 lines of code.
    return INFINITY; // currently reports no intersection
    // ----------- STUDENT CODE END ------------
}

float findIntersectionWithCone( Ray ray, vec3 center, vec3 apex, float radius, out Intersection out_intersect ) {
    vec3 axis   = center - apex;
    float len   = length( axis );
    axis = normalize( axis );

    // -- infinite cone
    Intersection intersect;
    float best_dist = INFINITY;
    float dist;

    // -- infinite cone
    dist = getIntersectOpenCone( ray, apex, axis, len, radius, intersect );
    chooseCloserIntersection( dist, best_dist, intersect, out_intersect );

    // -- caps
    dist = getIntersectDisc( ray, center, axis, radius, intersect );
    chooseCloserIntersection( dist, best_dist, intersect, out_intersect );

    return best_dist;
}

#define MAX_RECURSION 8

// Code copied directly from http://www.neilmendoza.com/glsl-rotation-about-an-arbitrary-axis/ as told to do so by Riley
mat4 rotationMatrix(vec3 axis, float angle)
{
    axis = normalize(axis);
    float s = sin(angle);
    float c = cos(angle);
    float oc = 1.0 - c;

    return mat4(oc * axis.x * axis.x + c,           oc * axis.x * axis.y - axis.z * s,  oc * axis.z * axis.x + axis.y * s,  0.0,
                oc * axis.x * axis.y + axis.z * s,  oc * axis.y * axis.y + c,           oc * axis.y * axis.z - axis.x * s,  0.0,
                oc * axis.z * axis.x - axis.y * s,  oc * axis.y * axis.z + axis.x * s,  oc * axis.z * axis.z + c,           0.0,
                0.0,                                0.0,                                0.0,                                1.0);
}

vec3 calculateSpecialDiffuseColor( Material mat, vec3 posIntersection, vec3 normalVector ) {
    // ----------- STUDENT CODE BEGIN ------------
    if ( mat.special == CHECKERBOARD ) {
        // do something here for checkerboard
        vec3 zAxis = vec3(0.0, 0.0, 1.0);
        vec3 normalizedAxis = cross(normalVector, zAxis);
        float angle = acos( dot(normalizedAxis, normalVector) );
        mat4 rotationMatrix = rotationMatrix(normalizedAxis, angle);

        vec4 oldCoordinates = vec4(posIntersection, 0.0);
        vec4 newCoordinates = oldCoordinates * rotationMatrix;
        float x = floor(newCoordinates.x);
        float y = floor(newCoordinates.y);
        float total = x + y;
        bool isEven = mod(total, 2.0) < EPS;

        vec3 blackColor = 0.5 * mat.color;
        vec3 whiteColor = 1.0 * mat.color;

        if (isEven) { return blackColor; }
        else        { return whiteColor; }
        // ----------- Our reference solution uses 21 lines of code.
    }
    else if ( mat.special == MYSPECIAL ) {
        // do something here for myspecial
        // ----------- Our reference solution uses 2 lines of code.
    }

    return mat.color; // special materials not implemented. just return material color.
    // ----------- STUDENT CODE END ------------
}

vec3 calculateDiffuseColor( Material mat, vec3 posIntersection, vec3 normalVector ) {
    // Special colors
    if ( mat.special != NONE ) {
        return calculateSpecialDiffuseColor( mat, posIntersection, normalVector );
    }
    return vec3( mat.color );
}

// check if position pos in in shadow with respect to a particular light.
// lightVec is the vector from that position to that light
bool pointInShadow( vec3 pos, vec3 lightVec ) {
    // ----------- STUDENT CODE BEGIN ------------
    Ray ray;
    ray.origin    = pos;
    ray.direction = normalize( lightVec );

    Material out_mat;
    Intersection out_intersect;
    float hitLength = rayIntersectScene( ray, out_mat, out_intersect );

    float distanceToIntersection = length(out_intersect.position - pos);
    // length > 0

    if (hitLength > -EPS && hitLength < length(lightVec) - EPS) { return true; }
    // ----------- Our reference solution uses 10 lines of code.
    return false;
    // ----------- STUDENT CODE END ------------
}

vec3 getLightContribution( Light light, Material mat, vec3 posIntersection, vec3 normalVector, vec3 eyeVector, bool phongOnly, vec3 diffuseColor ) {

    vec3 lightVector = light.position - posIntersection;

    if ( pointInShadow( posIntersection, lightVector ) ) {
        return vec3( 0.0, 0.0, 0.0 );
    }

    if ( mat.materialType == PHONGMATERIAL || mat.materialType == LAMBERTMATERIAL ) {
        vec3 contribution = vec3( 0.0, 0.0, 0.0 );

        // get light attenuation
        float dist = length( lightVector );
        float attenuation = light.attenuate * dist * dist;

        float diffuseIntensity = max( 0.0, dot( normalVector, lightVector ) ) * light.intensity;

        // glass and mirror objects have specular highlights but no diffuse lighting
        if ( !phongOnly ) {
            contribution += diffuseColor * diffuseIntensity * light.color / attenuation;
        }

        if ( mat.materialType == PHONGMATERIAL ) {
            // ----------- STUDENT CODE BEGIN ------------
            vec3 phongTerm = vec3( 0.0, 0.0, 0.0 ); // not implemented yet, so just add black
            // ----------- Our reference solution uses 10 lines of code.
            // ----------- STUDENT CODE END ------------
            contribution += phongTerm;
        }

        return contribution;
    }
    else {
        return diffuseColor;
    }

}

vec3 calculateColor( Material mat, vec3 posIntersection, vec3 normalVector, vec3 eyeVector, bool phongOnly ) {
	vec3 diffuseColor = calculateDiffuseColor( mat, posIntersection, normalVector );

	vec3 outputColor = vec3( 0.0, 0.0, 0.0 ); // color defaults to black when there are no lights

    for ( int i=0; i<MAX_LIGHTS; i++ ) {

        if( i>=numLights ) break; // because GLSL will not allow looping to numLights

        outputColor += getLightContribution( lights[i], mat, posIntersection, normalVector, eyeVector, phongOnly, diffuseColor );
	}

	return outputColor;
}

// find reflection or refraction direction ( depending on material type )
vec3 calcReflectionVector( Material material, vec3 direction, vec3 normalVector, bool isInsideObj ) {
    if( material.materialReflectType == MIRRORREFLECT ) {
        return reflect( direction, normalVector );
    }
    // the material is not mirror, so it's glass.
    // compute the refraction direction...

    // ----------- STUDENT CODE BEGIN ------------
    // see lecture 13 slide ( lighting ) on Snell's law
    // the eta below is eta_i/eta_r
    float eta = ( isInsideObj ) ? 1.0/material.refractionRatio : material.refractionRatio;
    // ----------- Our reference solution uses 11 lines of code.

    return reflect( direction, normalVector ); // return mirror direction so you can see something
    // ----------- STUDENT CODE END ------------
}

vec3 traceRay( Ray ray ) {
    Material hitMaterial;
    Intersection intersect;

    vec3 resColor  = vec3( 0.0, 0.0, 0.0 );
    vec3 resWeight = vec3( 1.0, 1.0, 1.0 );

    bool isInsideObj = false;

    for ( int depth = 0; depth < MAX_RECURSION; depth++ ) {

        float hit_length = rayIntersectScene( ray, hitMaterial, intersect );

        if ( hit_length < EPS || hit_length >= INFINITY ) break;

        vec3 posIntersection = intersect.position;
        vec3 normalVector    = intersect.normal;

        vec3 eyeVector = normalize( ray.origin - posIntersection );
        if ( dot( eyeVector, normalVector ) < 0.0 )
            { normalVector = -normalVector; isInsideObj = true; }
        else isInsideObj = false;

        bool reflective = ( hitMaterial.materialReflectType == MIRRORREFLECT ||
                            hitMaterial.materialReflectType == GLASSREFLECT );
		vec3 outputColor = calculateColor( hitMaterial, posIntersection, normalVector, eyeVector, reflective );

        float reflectivity = hitMaterial.reflectivity;

        // check to see if material is reflective ( or refractive )
        if ( !reflective || reflectivity < EPS ) {
            resColor += resWeight * outputColor;
            break;
        }

        // bounce the ray
        vec3 reflectionVector = calcReflectionVector( hitMaterial, ray.direction, normalVector, isInsideObj );
        ray.origin = posIntersection;
        ray.direction = normalize( reflectionVector );

        // add in the color of the bounced ray
        resColor += resWeight * outputColor;
        resWeight *= reflectivity;
    }

    return resColor;
}

void main( ) {
    float cameraFOV = 0.8;
    vec3 direction = vec3( v_position.x * cameraFOV * width/height, v_position.y * cameraFOV, 1.0 );

    Ray ray;
	  ray.origin    = vec3( uMVMatrix * vec4( camera, 1.0 ) );
    ray.direction = normalize( vec3( uMVMatrix * vec4( direction, 0.0 ) ) );

    // trace the ray for this pixel
    vec3 res = traceRay( ray );

    // paint the resulting color into this pixel
    gl_FragColor = vec4( res.x, res.y, res.z, 1.0 );
}
