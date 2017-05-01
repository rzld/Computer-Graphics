//RASTERISATION

//Progress 11/04/2017:
// lighting still not work!

#include <iostream>
#include <glm/glm.hpp>
#include <glm/fwd.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModel.h"
#include <limits>

#define PI 3.14159

using namespace std;
using glm::vec3;
using glm::mat3;
using glm::vec2;

/* ----------------------------------------------------------------------------*/
/* GLOBAL VARIABLES                                                            */

const int SCREEN_WIDTH = 500;
const int SCREEN_HEIGHT = 500;
SDL_Surface* screen;
int t;

vector<Triangle> triangles;
vec3 cameraPos(0, 0, -2.0);
float focalLength = 250.0;
mat3 R;
float yaw = 0;
float depthBuffer[SCREEN_HEIGHT][SCREEN_WIDTH];
vec3 currentColor;
vec3 lightPos(0, -0.5, -0.7);
vec3 lightPower = 1.1f*vec3(1,1,1);
vec3 indirectLight = 0.5f*vec3(1,1,1);
vec3 surfacePoint;

struct Pixel
{
	int x;
	int y;
	float zinv;
	vec3 illumination;
};
struct Vertex
{
	vec3 position;
	vec3 normal;
	vec3 reflectance;
};

/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

void Update();
void Draw();
void VertexShader(const vec3& v, glm::ivec2& p);
void Interpolate(glm::ivec2 a, glm::ivec2 b, vector<glm::ivec2>& result);
void DrawLineSDL(SDL_Surface* surface, glm::ivec2 a, glm::ivec2 b, vec3 color);
void DrawPolygonEdges(const vector<vec3>& vertices);
void ComputePolygonRows(const vector<glm::ivec2>& vertexPixels,
			      vector<glm::ivec2>& leftPixels,
			      vector<glm::ivec2>& rightPixels);
void DrawPolygonRows(const vector<glm::ivec2>& leftPixels, const vector<glm::ivec2>& rightPixels);
void DrawPolygon(const vector<vec3>& vertices);
//void DrawPolygonP(const vector<vec3>& vertices);
void DrawPolygonP(const vector<Vertex>& vertices);
void InterpolateP(Pixel a, Pixel b, vector<Pixel>& result);
void ComputePolygonRowsP(const vector<Pixel>& vertexPixels, vector<Pixel>& leftPixels, vector<Pixel>& rightPixels);
void DrawPolygonRowsP(const vector<Pixel>& leftPixels, const vector<Pixel>& rightPixels);
//void VertexShaderP(const vec3& v, Pixel& p);
void VertexShaderP(const Vertex& v, Pixel& p);
void PixelShader(const Pixel& p);
//void DrawPolygonEdgesP(const vector<vec3>& vertices);
void DrawPolygonEdgesP(const vector<Vertex>& vertices);

int main( int argc, char* argv[] )
{
	screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT );
	t = SDL_GetTicks();	// Set start value for timer.

	LoadTestModel(triangles);
	cout << "Load test model" << endl;

	// ComputePolygonRows test
	/*
	vector<glm::ivec2> vertexPixels(3);
	vertexPixels[0] = glm::ivec2(10, 5);
	vertexPixels[1] = glm::ivec2(5, 10);
	vertexPixels[2] = glm::ivec2(15, 15);
	vector<glm::ivec2> leftPixels;
	vector<glm::ivec2> rightPixels;
	ComputePolygonRows(vertexPixels, leftPixels, rightPixels);
	for (int row=0; row<leftPixels.size(); ++row)
	{
		cout << "Start: ("
		     << leftPixels[row].x << ","
		     << leftPixels[row].y << "). "
		     << "End: ("
		     << rightPixels[row].x << ","
		     << rightPixels[row].y << "). " << endl;
	}
	*/

	while( NoQuitMessageSDL() )
	{
		Update();
		Draw();
	}

	SDL_SaveBMP( screen, "screenshot.bmp" );
	return 0;
}

void Update()
{
	// Compute frame time:
	int t2 = SDL_GetTicks();
	float dt = float(t2-t);
	t = t2;
	cout << "Render time: " << dt << " ms." << endl;
}

void Draw()
{
	SDL_FillRect(screen,0,0);

	if( SDL_MUSTLOCK(screen) )
		SDL_LockSurface(screen);

	for (int y=0; y<SCREEN_HEIGHT; ++y)
		for (int x=0; x<SCREEN_WIDTH; ++x)
			depthBuffer[y][x] = 0;

	for (int i=0; i<(int)triangles.size(); ++i)
	{
		//cout << i << "/" << triangles.size() << endl;
		/*vector<vec3> vertices(3);
		vertices[0] = triangles[i].v0;
		vertices[1] = triangles[i].v1;
		vertices[2] = triangles[i].v2;*/

		vector<Vertex> vertices(3);
		vertices[0].position = vec3(triangles[i].v0);
		vertices[1].position = vec3(triangles[i].v1);
		vertices[2].position = vec3(triangles[i].v2);

		vertices[0].normal = vec3(triangles[i].normal);
		vertices[1].normal = vec3(triangles[i].normal);
		vertices[2].normal = vec3(triangles[i].normal);

		vertices[0].reflectance = vec3(triangles[i].color);
		vertices[1].reflectance = vec3(triangles[i].color);
		vertices[2].reflectance = vec3(triangles[i].color);

		surfacePoint = vec3(triangles[i].v0);

		currentColor = vec3(triangles[i].color);

		// dots only
		for (int v=0; v<3; ++v)
		{
			glm::ivec2 projPos;
			//VertexShader(vertices[v], projPos);
			VertexShader(vertices[v].position, projPos);
			//cout << projPos.x << " " << projPos.y << endl;
			//vec3 color(1.0, 1.0, 1.0);
			//PutPixelSDL(screen, projPos.x, projPos.y, color);
		}

		//DrawPolygonEdgesP(vertices);
		DrawPolygonP(vertices);
	}

	//original
	/*for( int y=0; y<SCREEN_HEIGHT; ++y )
	{
		for( int x=0; x<SCREEN_WIDTH; ++x )
		{
			vec3 color( 1.0, 0.0, 0.0 );
			PutPixelSDL( screen, x, y, color );
		}
	}*/

	if( SDL_MUSTLOCK(screen) )
		SDL_UnlockSurface(screen);

	SDL_UpdateRect( screen, 0, 0, 0, 0 );
}
void VertexShader(const vec3& v, glm::ivec2& p)
{
	// take 3D position of a vertex v and compute its 2D img position and store it in p
	vec3 pos = (v - cameraPos) * R;

	p.x = (focalLength * ( pos.x /  pos.z )) + (SCREEN_WIDTH / 2);
	p.y = (focalLength * ( pos.y /  pos.z )) + (SCREEN_HEIGHT / 2);
}

void Interpolate(glm::ivec2 a, glm::ivec2 b, vector<glm::ivec2>& result)
{
	int N = result.size();
	vec2 step = vec2(b-a) / float(max(N-1,1));
	vec2 current(a);

	for (int i=0; i<N; ++i)
	{
		result[i] = current;
		current += step;
	}
}

//draw line
void DrawLineSDL(SDL_Surface* surface, glm::ivec2 a, glm::ivec2 b, vec3 color)
{
	glm::ivec2 delta = glm::abs(a-b);
	int pixels = glm::max(delta.x, delta.y) + 1;

	vector<glm::ivec2> line(pixels);
	Interpolate(a, b, line);

	for (int i=0; i<(int)line.size(); i++)
	{
		PutPixelSDL(surface, line[i].x, line[i].y, color);
	}
}

//draw edges of triangle
void DrawPolygonEdges(const vector<vec3>& vertices)
{
	int V = vertices.size();

	//Transform each vertex from 3D world position to 2D image position
	vector<glm::ivec2> projectedVertices(V);
	for (int i=0; i<V; i++)
	{
		VertexShader(vertices[i], projectedVertices[i]);
	}

	//Loop over all vertices and draw the edge from it to the next vertex
	for (int i=0; i<V; i++)
	{
		int j = (i+1)%V; //next vertex
		vec3 color(1, 1, 1);
		DrawLineSDL(screen, projectedVertices[i], projectedVertices[j], color);
	}
}

//input: projPos
void ComputePolygonRows(const vector<glm::ivec2>& vertexPixels,
			      vector<glm::ivec2>& leftPixels,
			      vector<glm::ivec2>& rightPixels)
{
	int minY = 999;
	int maxY = -999;
	int rows;

	// 1. Find max & min y-value of the polygon and compute the number of rows it occupies
	//cout << "Step 1" << endl;
	for (int i=0; i<(int)vertexPixels.size(); i++)
	{
		minY = glm::min(minY, vertexPixels[i].y);
		maxY = glm::max(maxY, vertexPixels[i].y);
	}
	rows = maxY - minY + 1;
	//cout << maxY << "-" << minY << "=" << rows << endl;

	// 2. Resize leftPixels and rightPixels so that they have an element for each row
	//cout << "Step 2, " << endl;
	leftPixels.resize(rows);
	rightPixels.resize(rows);

	// 3. Initialize the x-coordinates in leftPixels to some really large value
	//    and the x-coordinates in the rightPixels to some really small values
	//cout << "Step 3" << endl;
	for (int i=0; i<rows; i++)
	{
		leftPixels[i].x = +numeric_limits<int>::max();
		rightPixels[i].x = -numeric_limits<int>::max();

		leftPixels[i].y = minY + i;
		rightPixels[i].y = minY + i;
	}

	// 4. Loop through all edges of the polygon and use linear interpolation to find
	//    the x-coordinate for each row it occupies. Update the corresponding values in rightPixels
	//    and leftPixels
	//cout << "Step 4" << endl;
	for (int i=0; i<(int)vertexPixels.size(); i++)
	{
		int j = (i+1)%vertexPixels.size(); //next vertex
		glm::ivec2 delta = glm::abs(vertexPixels[i]-vertexPixels[j]);
		int pixels = glm::max(delta.x, delta.y) + 1;

		vector<glm::ivec2> line(pixels);
		Interpolate(vertexPixels[i], vertexPixels[j], line);

		// check for every point in line
		for (int a=0; a<(int)line.size(); a++)
		{
			//cout << line[i].x << " " << line[i].y << endl;
			// check for every point in left & right pixels
			for (int b=0; b<rows; b++)
			{
				// if they are in the same y line
				if (leftPixels[b].y == line[a].y)
				{
					int leftX = glm::min(leftPixels[b].x, line[a].x);
					int rightX = glm::max(rightPixels[b].x, line[a].x);

					leftPixels[b].x = leftX;
					rightPixels[b].x = rightX;
				}
			}
		}
	}
}

void DrawPolygonRows(const vector<glm::ivec2>& leftPixels, const vector<glm::ivec2>& rightPixels)
{
	for (int i=0; i<(int)leftPixels.size(); i++)
	{
		vector<glm::ivec2> rowPixels((rightPixels[i].x - leftPixels[i].x)+1);
		Interpolate(leftPixels[i], rightPixels[i], rowPixels);
		for (int i=0; i<(int)rowPixels.size(); i++)
		{
			PutPixelSDL(screen, rowPixels[i].x, rowPixels[i].y, currentColor);
		}
	}
}

void DrawPolygon(const vector<vec3>& vertices)
{
	int V = vertices.size();

	vector<glm::ivec2> vertexPixels(V);
	for (int i=0; i<V; i++)
	{
		VertexShader(vertices[i], vertexPixels[i]);
	}

	vector<glm::ivec2> leftPixels;
	vector<glm::ivec2> rightPixels;
	ComputePolygonRows(vertexPixels, leftPixels, rightPixels);
	for (int row=0; row<(int)leftPixels.size(); ++row)
	{
		cout << "Start: ("
		     << leftPixels[row].x << ","
		     << leftPixels[row].y << "). "
		     << "End: ("
		     << rightPixels[row].x << ","
		     << rightPixels[row].y << "). " << endl;
	}
	DrawPolygonRows(leftPixels, rightPixels);
}

//void VertexShaderP(const vec3& v, Pixel& p)
void VertexShaderP(const Vertex& v, Pixel& p)
{
	// take 3D position of a vertex v and compute its 2D img position and store it in p
	vec3 pos = (v.position - cameraPos) * R;

	p.zinv = 1/pos.z;
	p.x = (focalLength * (pos.x / pos.z)) + (SCREEN_WIDTH / 2);
	p.y = (focalLength * (pos.y / pos.z)) + (SCREEN_HEIGHT / 2);
	//cout << p.x << " " << p.y << " " << p.zinv << endl;
	//depthBuffer[p.y][p.x] = p.zinv;

	//light
	float radius;
	float area;
	vec3 _n, _r, d;

	_n = v.normal;
	radius = glm::distance(lightPos, surfacePoint);
	area = 4 * PI * radius * radius;

	//cout << radius << " " << area << endl;

	_r = glm::normalize(lightPos - surfacePoint);

	float pMax = glm::max((float)glm::dot(_r, _n), (float)0.0);

	//cout << pMax << endl;

	if (pMax > 0.0)
	{
		d = lightPower/area;
	}
	else
	{
		d = vec3(0.0, 0.0, 0.0);
	}

	//d = lightPower/area;
	p.illumination = v.reflectance * (d + indirectLight);

	//cout << p.illumination.x << " " << p.illumination.y << " " << p.illumination.z << endl;
}

//void DrawPolygonP(const vector<vec3>& vertices)
void DrawPolygonP(const vector<Vertex>& vertices)
{
	int V = vertices.size();

	vector<Pixel> vertexPixels(V);
	for (int i=0; i<V; i++)
	{
		VertexShaderP(vertices[i], vertexPixels[i]);
	}

	vector<Pixel> leftPixels;
	vector<Pixel> rightPixels;
	ComputePolygonRowsP(vertexPixels, leftPixels, rightPixels);
	DrawPolygonRowsP(leftPixels, rightPixels);
}

void InterpolateP(Pixel a, Pixel b, vector<Pixel>& result)
{
	int N = result.size();

	float stepX = (b.x - a.x) / float(glm::max(N-1, 1));
	float stepY = (b.y - a.y) / float(glm::max(N-1, 1));
	float stepZ = (b.zinv - a.zinv) / float(glm::max(N-1, 1));
	vec3 stepI = (b.illumination - a.illumination) / float(glm::max(N-1, 1));

	//cout << a.illumination.x << " " << a.illumination.y << " " << a.illumination.z << endl;
	//cout << b.illumination.x << " " << b.illumination.y << " " << b.illumination.z << endl << endl;
	//cout << stepI.x << " " << stepI.y << " " << stepI.z << endl;

	float currentX = (float)a.x;
	float currentY = (float)a.y;
	float currentZ = a.zinv;
	vec3 currentI = a.illumination;

	for (int i=0; i<N; i++)
	{
		result[i].x = (int)currentX;
		result[i].y = (int)currentY;
		result[i].zinv = currentZ;
		result[i].illumination = currentI;

		currentX += stepX;
		currentY += stepY;
		currentZ += stepZ;
		currentI += stepI;
	}
}

void ComputePolygonRowsP(const vector<Pixel>& vertexPixels,
			      vector<Pixel>& leftPixels,
			      vector<Pixel>& rightPixels)
{
	int minY = 999;
	int maxY = -999;
	int rows;

	// 1. Find max & min y-value of the polygon and compute the number of rows it occupies
	for (int i=0; i<(int)vertexPixels.size(); i++)
	{
		minY = glm::min(minY, vertexPixels[i].y);
		maxY = glm::max(maxY, vertexPixels[i].y);
	}
	rows = maxY - minY + 1;

	// 2. Resize leftPixels and rightPixels so that they have an element for each row
	leftPixels.resize(rows);
	rightPixels.resize(rows);

	// 3. Initialize the x-coordinates in leftPixels to some really large value
	//    and the x-coordinates in the rightPixels to some really small values
	for (int i=0; i<rows; i++)
	{
		leftPixels[i].x = +numeric_limits<int>::max();
		rightPixels[i].x = -numeric_limits<int>::max();

		leftPixels[i].y = minY + i;
		rightPixels[i].y = minY + i;
	}

	// 4. Loop through all edges of the polygon and use linear interpolation to find
	//    the x-coordinate for each row it occupies. Update the corresponding values in rightPixels
	//    and leftPixels
	for (int i=0; i<(int)vertexPixels.size(); i++)
	{
		int j = (i+1) % vertexPixels.size(); //next vertex
		//Pixel delta = glm::abs(vertexPixels[i]-vertexPixels[j]);
		int deltaX = glm::abs(vertexPixels[i].x - vertexPixels[j].x);
		int deltaY = glm::abs(vertexPixels[i].y - vertexPixels[j].y);
		int pixels = glm::max(deltaX, deltaY) + 1;

		vector<Pixel> line;
		line.resize(pixels);
		InterpolateP(vertexPixels[i], vertexPixels[j], line);

		// check for every point in line
		for (int a=0; a<(int)line.size(); a++)
		{
			//cout << line[a].x << " " << line[a].y << " " << line[a].zinv << endl;
			// check for every point in left & right pixels
			for (int b=0; b<rows; b++)
			{
				// if they are in the same y line
				//how about the zinv?
				if (leftPixels[b].y == line[a].y)
				{
					int leftX = glm::min(leftPixels[b].x, line[a].x);
					int rightX = glm::max(rightPixels[b].x, line[a].x);

					//leftPixels[b].x = leftX;
					//rightPixels[b].x = rightX;

					if (line[a].x < leftPixels[b].x)
					{
						leftPixels[b].x = leftX;
						leftPixels[b].zinv = line[a].zinv;
						leftPixels[b].illumination = vec3(line[a].illumination);
						//depthBuffer[leftPixels[b].y][leftPixels[b].x] = leftPixels[b].zinv;
					}
					if (line[a].x > rightPixels[b].x)
					{
						rightPixels[b].x = rightX;
						rightPixels[b].zinv = line[a].zinv;
						rightPixels[b].illumination = vec3(line[a].illumination);
						//depthBuffer[rightPixels[b].y][rightPixels[b].x] = rightPixels[b].zinv;
					}
				}
			}
		}
	}
}

void DrawPolygonRowsP(const vector<Pixel>& leftPixels, const vector<Pixel>& rightPixels)
{
	for (int i=0; i<(int)leftPixels.size(); i++)
	{
		vector<Pixel> rowPixels((rightPixels[i].x - leftPixels[i].x)+1);
		InterpolateP(leftPixels[i], rightPixels[i], rowPixels);

		for (int j=0; j<(int)rowPixels.size(); j++)
		{
			PixelShader(rowPixels[j]);
			//cout << rowPixels[j].illumination.x << " " << rowPixels[j].illumination.y << " " << rowPixels[j].illumination.z << endl;
		}
	}
}

//void DrawPolygonEdgesP(const vector<vec3>& vertices)
void DrawPolygonEdgesP(const vector<Vertex>& vertices)
{
	int V = vertices.size();

	//Transform each vertex from 3D world position to 2D image position
	vector<Pixel> projectedVertices(V);
	for (int i=0; i<V; i++)
	{
		VertexShaderP(vertices[i], projectedVertices[i]);
	}

	//Loop over all vertices and draw the edge from it to the next vertex
	for (int i=0; i<V; i++)
	{
		int j = (i+1)%V; //next vertex
		glm::ivec2 startV(projectedVertices[i].x, projectedVertices[i].y);
		glm::ivec2 endV(projectedVertices[j].x, projectedVertices[j].y);
		vec3 color(1, 1, 1);
		DrawLineSDL(screen, startV, endV, currentColor);
	}
}

void PixelShader(const Pixel& p)
{
	//cout << p.illumination.x << " " << p.illumination.y << " " << p.illumination.z << endl;

	if(p.zinv > depthBuffer[p.y][p.x])
	{
		depthBuffer[p.y][p.x] = p.zinv;
		PutPixelSDL(screen, p.x, p.y, p.illumination);
	}
}
