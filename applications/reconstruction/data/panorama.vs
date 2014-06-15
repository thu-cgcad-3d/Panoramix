
uniform mat4 worldViewProj;


varying vec3 pixelPosition;
varying vec3 pixelNormal;

void main(void)
{
    pixelPosition = gl_Vertex.xyz / gl_Vertex.w;
    pixelNormal = gl_Normal;
    gl_Position = worldViewProj * gl_Vertex;
}