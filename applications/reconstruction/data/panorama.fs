
uniform sampler2D Image;
uniform vec3 CasterCenter;

varying vec3 pixelPosition;
varying vec3 pixelNormal;

void main(void)
{
    vec3 direction = pixelPosition - CasterCenter;
    float longi = atan(direction.z, direction.x);
    float lati = - asin(direction.y / length(direction));
    vec2 texCoord = vec2(longi / 3.1415926535897932 / 2.0 + 0.5, lati / 3.1415926535897932 + 0.5);
    gl_FragColor = texture2D(Image, texCoord);
   	//gl_FragColor = vec4(1.0, 0.0, 0.0, 1.0);
}