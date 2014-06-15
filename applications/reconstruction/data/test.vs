//#version 150

uniform mat4 worldViewProj;

//in vec4 vertex;

//out vec4 position;

void main(void)
{
    //position = worldViewProj * vertex;
    gl_Position = worldViewProj * gl_Vertex;
}