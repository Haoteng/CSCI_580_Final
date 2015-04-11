/* CS580 Homework 3 */

#include	<math.h>
#include	"gz.h"
#include	"rend.h"
#include	"PooyanToolset.h"
#include	"disp.h"

void rasterize(GzRender *, GzCoord, GzCoord, GzTextureIndex, GzCoord, GzCoord,
		GzTextureIndex, GzCoord, GzCoord, GzTextureIndex);
void calColor(GzRender *, GzCoord, GzColor &, GzColor);

int GzNewRender(GzRender **render, GzRenderClass renderClass,
		GzDisplay *display) {
	/*
	 - malloc a renderer struct
	 - keep closed until all inits are done
	 - setup Xsp and anything only done once
	 - save pointer to display for init and rasterizer
	 - check for legal class GZ_Z_BUFFER_RENDER
	 - init default camera - normalize "up vector" just in case.
	 */
	if (!display)
		return GZ_FAILURE;
	if (renderClass != GZ_Z_BUFFER_RENDER) // - check for legal class GZ_Z_BUFFER_RENDER
		return GZ_FAILURE;
	(*render) = new GzRender(); //- malloc a renderer struct
	(*render)->renderClass = renderClass;
	(*render)->matlevel = -1; // empty stack
	double d = 1 / tan(PooyanToolset::degreeToRadian(DEFAULT_FOV / 2)); // 1/d=tan(FOV/2)
	/*
	 Xsp = xs/2 0 0 xs/2
	 0 -ys/2 0 ys/2
	 0 0 Zmax/d 0
	 0 0 0 1
	 */
	GzMatrix Xsp = { { display->xres / 2, 0, 0, display->xres / 2 }, { 0,
			-display->yres / 2, 0, display->yres / 2 },
			{ 0, 0, MAX_INT / d, 0 }, { 0, 0, 0, 1 } };
	PooyanToolset::copyFromVoid((*render)->Xsp, Xsp); // - setup Xsp and anything only done once
	GzCamera defaultCamera;
	defaultCamera.FOV = DEFAULT_FOV;
	defaultCamera.lookat[X] = 0; /* default look-at point = 0,0,0 */
	defaultCamera.lookat[Y] = 0; /* default look-at point = 0,0,0 */
	defaultCamera.lookat[Z] = 0; /* default look-at point = 0,0,0 */
	defaultCamera.position[X] = DEFAULT_IM_X;
	defaultCamera.position[Y] = DEFAULT_IM_Y;
	defaultCamera.position[Z] = DEFAULT_IM_Z;
	defaultCamera.worldup[X] = 0;
	defaultCamera.worldup[Y] = 1;
	defaultCamera.worldup[Z] = 0;
	(*render)->camera = defaultCamera; //init default camera - normalize "up vector" just in case.
	(*render)->display = display; //- save pointer to display for init and rasterizer
	(*render)->open = 1; // - keep closed until all inits are done
	(*render)->numlights = 0; // no light to begin with
	(*render)->shift[X]=0;
	(*render)->shift[Y]=0;
	return GZ_SUCCESS;
}

int GzFreeRender(GzRender *render) {
	/*
	 -free all renderer resources
	 -it doesn't exist anymore
	 */

	if (!render)
		return GZ_FAILURE;
	delete render; //-free all renderer resources
	return GZ_SUCCESS;
}

int GzBeginRender(GzRender *render) {
	/*
	 - set up for start of each frame - clear frame buffer
	 - compute Xiw and projection xform Xpi from camera definition
	 - init Xform stack - put Xsp at base of stack, push on Xpi and Xiw
	 - now stack contains Xsw and app can push model Xforms when it want to
	 */
	int index[] = { X, Y, Z };
	if (!render)
		return GZ_FAILURE;
	GzInitDisplay(render->display); //- set up for start of each frame - clear frame buffer
	GzCamera camera = render->camera;
	double d = 1 / tan(PooyanToolset::degreeToRadian(camera.FOV / 2)); // 1/d=tan(FOV/2)
	render->Xsp[2][2] = MAX_INT / d;
	/*
	 Xpi = 1 0 0 0
	 0 1 0 0
	 0 0 1 0
	 0 0 1/d 1
	 */
	// -- projection xform Xpi from camera definition
	GzMatrix Xpi = { { 1, 0, 0, 0 }, { 0, 1, 0, 0 }, { 0, 0, 1, 0 }, { 0, 0, 1
			/ d, 1 } };
	GzCoord z;
	for (int i = 0; i < 3; i++) // this also can be done by operator overloading, do this next time
		z[index[i]] = camera.lookat[index[i]] - camera.position[index[i]];
	if (PooyanToolset::normalize(z)) // if the length of the vector is zero, it cannot be normalized
		return GZ_FAILURE;
	if (PooyanToolset::isEqual(camera.worldup, z)) // worldup and camera normalized vector cannot be equal
		return GZ_FAILURE;
	//  up' = up - (up * Z)Z
	float innerProductCameraWorldUp = PooyanToolset::innerProduct(z,
			camera.worldup);
	GzCoord y;
	for (int i = 0; i < 3; i++)
		y[index[i]] = camera.worldup[index[i]]
				- innerProductCameraWorldUp * z[index[i]];
	if (PooyanToolset::normalize(y))
		return GZ_FAILURE; // normalized worldUpP is y in fact s
	GzCoord x;
	PooyanToolset::crossProduct(y, z, x); // X = (Y x Z)
	//	- compute Xiw
	GzMatrix Xiw = { { x[X], x[Y], x[Z], -PooyanToolset::innerProduct(x,
			camera.position) }, { y[X], y[Y], y[Z],
			-PooyanToolset::innerProduct(y, camera.position) }, { z[X], z[Y],
			z[Z], -PooyanToolset::innerProduct(z, camera.position) }, { 0, 0, 0,
			1 } };
	// - init Xform stack - put Xsp at base of stack, push on Xpi and Xiw
	GzPushMatrix(render, render->Xsp);
	GzPushMatrix(render, Xpi);
	GzPushMatrix(render, Xiw);
	//- now stack contains Xsw and app can push model Xforms when it want to
	return GZ_SUCCESS;
}

int GzPutAttribute(GzRender *render, int numAttributes, GzToken *nameList,
		GzPointer *valueList) /* void** valuelist */
		{
	/*
	 - set renderer attribute states (e.g.: GZ_RGB_COLOR default color)
	 - later set shaders, interpolaters, texture maps, and lights
	 */
	if (!render)
		return GZ_FAILURE;
	for (int i = 0; i < numAttributes; i++)
		switch (nameList[i]) {
		case GZ_RGB_COLOR:
			PooyanToolset::copyFromVoid(render->flatcolor, valueList[i]);
			break;
		case GZ_INTERPOLATE:
			PooyanToolset::copyFromVoid(render->interp_mode, valueList[i]);
			break;
		case GZ_DIRECTIONAL_LIGHT:
			PooyanToolset::copyFromVoid(render->lights[render->numlights],
					valueList[i]);
			render->numlights++;
			break;
		case GZ_AMBIENT_LIGHT:
			PooyanToolset::copyFromVoid(render->ambientlight, valueList[i]);
			break;
		case GZ_AMBIENT_COEFFICIENT:
			PooyanToolset::copyFromVoid(render->Ka, valueList[i]);
			break;
		case GZ_DIFFUSE_COEFFICIENT:
			PooyanToolset::copyFromVoid(render->Kd, valueList[i]);
			break;
		case GZ_SPECULAR_COEFFICIENT:
			PooyanToolset::copyFromVoid(render->Ks, valueList[i]);
			break;
		case GZ_DISTRIBUTION_COEFFICIENT:
			PooyanToolset::copyFromVoid(render->spec, valueList[i]);
			break;
		case GZ_TEXTURE_MAP:
			render->tex_fun = (int (*)(float, float, GzColor))valueList[i];
			break;
		case GZ_AASHIFTX:
			PooyanToolset::copyFromVoid(render->shift[X], valueList[i]) ;
			break ;
		case GZ_AASHIFTY:
			PooyanToolset::copyFromVoid(render->shift[Y], valueList[i]) ;
			break ;
		default:
			return GZ_FAILURE;
			break;
		}

	return GZ_SUCCESS;

}

int GzPutTriangle(GzRender *render, int numParts, GzToken *nameList,
		GzPointer *valueList)
		/* numParts - how many names and values */
		{
	/*
	 - pass in a triangle description with tokens and values corresponding to
	 GZ_NULL_TOKEN:		do nothing - no values
	 GZ_POSITION:		3 vert positions in model space
	 - Invoke the scan converter and return an error code
	 */
	float vertex[3][3];
	float normal[3][3];
	float puv[3][2];
	if (!render)
		return GZ_FAILURE;
	for (int i = 0; i < numParts; i++) {
		switch (nameList[i]) {
		case GZ_NULL_TOKEN:
			continue;
			break;
		case GZ_POSITION:
			PooyanToolset::copyFromVoid(vertex, valueList[i]);
			// - Xform vert coordinates
			for (int ii = 0; ii < 3; ii++) {
				float xp[4] = { vertex[ii][0], vertex[ii][1], vertex[ii][2], 1 };
				float xframed[4] = { 0, 0, 0, 0 };
				for (int jj = 0; jj < 4; jj++)
					for (int k = 0; k < 4; k++)
						xframed[jj] += render->Ximage[render->matlevel][jj][k]
								* xp[k];
				if (xframed[2] <= -1) //- Clip - discard any triangle with verts behind view plane (Z < 0)
					return GZ_FAILURE;
				for (int j = 0; j < 3; j++)
					vertex[ii][j] = xframed[j] / xframed[3];
			}
			break;
		case GZ_NORMAL:
			PooyanToolset::copyFromVoid(normal, valueList[i]);
			for (int nc = 0; nc < 3; nc++) { //normal
				float temp[4] =
						{ normal[nc][0], normal[nc][1], normal[nc][2], 1 };

				for (int row = 0; row < 3; row++) {
					normal[nc][row] = 0;

					for (int inner = 0; inner < 4; inner++) {
						normal[nc][row] +=
								render->Xnorm[(render->matlevel - 2)][row][inner]
										* temp[inner];

					}

				}

			}
			break;
		case GZ_TEXTURE_INDEX:
			PooyanToolset::copyFromVoid(puv, valueList[i]);
			break;
		default:
			return GZ_FAILURE;
		}
	}
	//Anti Alising
	for (int ii =0; ii<3; ii++ ){
		vertex[ii][X]+=render->shift[X];
		vertex[ii][Y]+=render->shift[Y];
	}
	//perspective correction
	for (int ii = 0; ii < 3; ii++) {
		float vz = vertex[ii][2] / (MAX_INT - vertex[ii][2]);
		puv[ii][0] = puv[ii][0] / (vz + 1);
		puv[ii][1] = puv[ii][1] / (vz + 1);
	}

	float y[] = { vertex[0][1], vertex[1][1], vertex[2][1] };

	int comb[][3] = { { 0, 1, 2 }, { 0, 2, 1 }, { 1, 0, 2 }, { 1, 2, 0 }, { 2,
			0, 1 }, { 2, 1, 0 } };

	for (int i = 0; i < 6; i++)
		if (y[comb[i][0]] < y[comb[i][1]] && y[comb[i][1]] < y[comb[i][2]]) {
			rasterize(render, vertex[comb[i][0]], normal[comb[i][0]],
					puv[comb[i][0]], vertex[comb[i][1]], normal[comb[i][1]],
					puv[comb[i][1]], vertex[comb[i][2]], normal[comb[i][2]],
					puv[comb[i][2]]);
			break;
		}

	return GZ_SUCCESS;
}

/* NOT part of API - just for general assistance */

short ctoi(float color) /* convert float color to GzIntensity short */
{
	return (short) ((int) (color * ((1 << 12) - 1)));
}

void rasterize(GzRender *render, GzCoord vMin, GzCoord nMin,
		GzTextureIndex uvMin, GzCoord vMid, GzCoord nMid, GzTextureIndex uvMid,
		GzCoord vMax, GzCoord nMax, GzTextureIndex uvMax) {

	float zero = 0;
	GzCoord minBoundry, maxBoundry = { 0 };
	minBoundry[0] = MIN(MIN(vMin[0],vMid[0]), vMax[0]);
	PooyanToolset::fitValue(minBoundry[0], zero, (float) render->display->xres);
	minBoundry[1] = vMin[1];
	PooyanToolset::fitValue(minBoundry[1], zero, (float) render->display->yres);
	maxBoundry[0] = MAX(MAX(vMin[0],vMid[0]), vMax[0]);
	PooyanToolset::fitValue(maxBoundry[0], zero, (float) render->display->xres);
	maxBoundry[1] = vMax[1];
	PooyanToolset::fitValue(maxBoundry[1], zero, (float) render->display->yres);

	float lineCoef[3][3] = { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } }; // coef[0]:min,mid coef[1]:mid,max coef[2]:max,min
	float vertex[3][3] = { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } };
	float normal[3][3] = { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } };
	float uv[3][2] = { { 0, 0 }, { 0, 0 }, { 0, 0 } };
	for (int i = 0; i < 3; i++) {
		vertex[0][i] = vMin[i];
		normal[0][i] = nMin[i];
		vertex[1][i] = vMid[i];
		normal[1][i] = nMid[i];
		vertex[2][i] = vMax[i];
		normal[2][i] = nMax[i];
	}
	for (int i = 0; i < 2; i++) {
		uv[0][i] = uvMin[i];
		uv[1][i] = uvMid[i];
		uv[2][i] = uvMax[i];
	}
	for (int i = 0; i < 3; i++) {
		PooyanToolset::lineEquation(lineCoef[i], vertex[i],
				vertex[(i + 1) % 3]);
	}
	float xIntersection = -1 * (lineCoef[2][1] * vMid[1] + lineCoef[2][2])
			/ lineCoef[2][0];
	int sign = SIGN(xIntersection - vMid[0]);

	float plateCoef[4] = { 0 };
	PooyanToolset::plateEquation(plateCoef, vertex[0], vertex[1], vertex[2]);
	float vCopy[3][3];
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			vCopy[i][j] = vertex[i][j];

	float normalCoef[3][4];
	float colorCoef[3][4];
	if (render->interp_mode == GZ_NORMALS)
		for (int i = 0; i < 3; i++) {
			for (int n = 0; n < 3; n++)
				vCopy[n][2] = normal[n][i];
			PooyanToolset::plateEquation(normalCoef[i], vCopy[0], vCopy[1],
					vCopy[2]);
		}
	else if (render->interp_mode == GZ_COLOR) {
		GzColor color[3];
		for (int i = 0; i < 3; i++) {
			GzColor texture;
			render->tex_fun(uv[i][0], uv[i][1], texture);
			calColor(render, normal[i], color[i], texture);
		}
		for (int i = 0; i < 3; i++) {
			for (int n = 0; n < 3; n++)
				vCopy[n][2] = color[n][i];
			PooyanToolset::plateEquation(colorCoef[i], vCopy[0], vCopy[1],
					vCopy[2]);
		}
	}
	float uvCoef[2][4];
	for (int i = 0; i < 2; i++) {
		for (int n = 0; n < 3; n++)
			vCopy[n][2] = uv[n][i];
		PooyanToolset::plateEquation(uvCoef[i], vCopy[0], vCopy[1], vCopy[2]);
	}
	for (int j = minBoundry[1]; j < maxBoundry[1]; j++)
		for (int i = minBoundry[0]; i < maxBoundry[0]; i++) {
			float lee[3] = { 0 };
			int signSum = 0;
			for (int ii = 0; ii < 3; ii++) {
				lee[ii] = lineCoef[ii][0] * i + lineCoef[ii][1] * j
						+ lineCoef[ii][2];
				signSum += sign * SIGN(lee[ii]);
			}
			if (signSum == -3) {
				GzIntensity dummy;
				GzDepth cameraZ;
				GzGetDisplay(render->display, i, j, &dummy, &dummy, &dummy,
						&dummy, &cameraZ);
				GzDepth z = PooyanToolset::calZ(plateCoef, i, j);
				if (cameraZ > z) { // z buffer checking
					// Calculating u and v considering i and j
					GzTextureIndex uvij;
					float vz = ((float) z) / (MAX_INT - z);
					for (int ii = 0; ii < 2; ii++) {
						uvij[ii] = PooyanToolset::calZ(uvCoef[ii], i, j);
						uvij[ii] *= (vz + 1); // image space
					}

					GzColor texture;

					render->tex_fun(uvij[X], uvij[Y], texture);

					// Color
					GzColor c;

					if (render->interp_mode == GZ_NORMALS) {
						GzCoord n;
						for (int ii = 0; ii < 3; ii++)
							n[ii] = PooyanToolset::calZ(normalCoef[ii], i, j);
						calColor(render, n, c, texture);
					} else if (render->interp_mode == GZ_COLOR)
						for (int ii = 0; ii < 3; ii++)
							c[ii] = PooyanToolset::calZ(colorCoef[ii], i, j)
									* texture[ii];
					else
						for (int ii = 0; ii < 3; i++)
							c[ii] = render->flatcolor[ii];

					GzIntensity giColor[3];
					for (int ii = 0; ii < 3; ii++) {
						PooyanToolset::fitValue(c[ii], (float) 0, (float) 1);
						giColor[ii] = ctoi(c[ii]);
					}

					GzPutDisplay(render->display, i, j, giColor[0], giColor[1],
							giColor[2], 1, z);

				}

			}
		}
}

//-------------------

int GzRotXMat(float degree, GzMatrix mat) {
// Create rotate matrix : rotate along x axis
// Pass back the matrix using mat value
	float r = PooyanToolset::degreeToRadian(degree);
	GzMatrix rx = { { 1, 0, 0, 0 }, { 0, cos(r), -sin(r), 0 }, { 0, sin(r), cos(
			r), 0 }, { 0, 0, 0, 1 } };
	PooyanToolset::copyFromVoid(mat, rx);
	return GZ_SUCCESS;
}

int GzRotYMat(float degree, GzMatrix mat) {
// Create rotate matrix : rotate along y axis
// Pass back the matrix using mat value
	float r = PooyanToolset::degreeToRadian(degree);
	GzMatrix ry = { { cos(r), 0, sin(r), 0 }, { 0, 1, 0, 0 }, { -sin(r), 0, cos(
			r), 0 }, { 0, 0, 0, 1 } };
	PooyanToolset::copyFromVoid(mat, ry);
	return GZ_SUCCESS;
}

int GzRotZMat(float degree, GzMatrix mat) {
// Create rotate matrix : rotate along z axis
// Pass back the matrix using mat value
	float r = PooyanToolset::degreeToRadian(degree);
	GzMatrix rx = { { cos(r), -sin(r), 0, 0 }, { sin(r), cos(r), 0, 0 }, { 0, 0,
			1, 0 }, { 0, 0, 0, 1 } };
	PooyanToolset::copyFromVoid(mat, rx);
	return GZ_SUCCESS;
}

int GzTrxMat(GzCoord translate, GzMatrix mat) {
// Create translation matrix
// Pass back the matrix using mat value
	GzMatrix t = { { 1, 0, 0, translate[X] }, { 0, 1, 0, translate[Y] }, { 0, 0,
			1, translate[Z] }, { 0, 0, 0, 1 } };
	PooyanToolset::copyFromVoid(mat, t);
	return GZ_SUCCESS;
}

int GzScaleMat(GzCoord scale, GzMatrix mat) {
// Create scaling matrix
// Pass back the matrix using mat value

	GzMatrix s = { { scale[X], 0, 0, 0 }, { 0, scale[Y], 0, 0 }, { 0, 0,
			scale[Z], 0 }, { 0, 0, 0, 1 } };
	PooyanToolset::copyFromVoid(mat, s);
	return GZ_SUCCESS;
}

int GzPutCamera(GzRender *render, GzCamera *camera) {
	/*
	 - overwrite renderer camera structure with new camera definition
	 - normalize "up vector" just in case (sanity check)
	 */
	if (!render)
		return GZ_FAILURE;
	PooyanToolset::copyFromVoid(render->camera, camera);
	double d = 1 / PooyanToolset::degreeToRadian(camera->FOV / 2);
	render->Xsp[2][2] = MAX_INT / d; // Zmax/d
	GzCamera& cam = render->camera;
	PooyanToolset::normalize(cam.worldup); // - normalize "up vector" just in case (sanity check)
	return GZ_SUCCESS;
}

int GzPushMatrix(GzRender *render, GzMatrix matrix) {
	/*
	 - push a matrix onto the Ximage stack
	 - check for stack overflow
	 */
	if (!render)
		return GZ_FAILURE;
	if (render->matlevel >= (MATLEVELS - 1))
		// stack overflow
		return GZ_FAILURE;
	if (render->matlevel < 0) {
		// an empty statck
		render->matlevel++;
		PooyanToolset::copyFromVoid(render->Ximage[render->matlevel], matrix);
	} else {
		render->matlevel++;
		// Ximage
		// Top of stack = Q1 = (Q0B) = (AB)
		for (int i = 0; i < 4; i++)
			for (int j = 0; j < 4; j++) {
				render->Ximage[render->matlevel][i][j] = 0;
				for (int k = 0; k < 4; k++)
					render->Ximage[render->matlevel][i][j] +=
							render->Ximage[render->matlevel - 1][i][k]
									* matrix[k][j];
			}
		// XNorm
		///Put zeros in right column (3 upper elements) prior to push
		GzMatrix copyMatrix;
		for (int i = 0; i < 4; i++)
			for (int j = 0; j < 4; j++)
				copyMatrix[i][j] = matrix[i][j];
		for (int i = 0; i < 3; i++)
			copyMatrix[i][3] = 0;
		// use any row/col and compute scale factor
		float K = 1 / PooyanToolset::lengthCalculation(matrix[1]);
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++)
				copyMatrix[i][j] *= K;
		}

		int xnormLevel = render->matlevel - 2;
		if (xnormLevel == 0) // the stack used to be empty
			PooyanToolset::copyFromVoid(render->Xnorm[xnormLevel], copyMatrix);
		else
			for (int i = 0; i < 4; i++)
				for (int j = 0; j < 4; j++) {
					render->Xnorm[xnormLevel][i][j] = 0;
					for (int k = 0; k < 4; k++)
						render->Xnorm[xnormLevel][i][j] +=
								render->Xnorm[xnormLevel - 1][i][k]
										* copyMatrix[k][j];
				}
	}

	return GZ_SUCCESS;
}

int GzPopMatrix(GzRender *render) {
	/*
	 - pop a matrix off the Ximage stack
	 - check for stack underflow
	 */
	if (!render)
		return GZ_FAILURE;
	if (render->matlevel < 0)
		//empty stack
		return GZ_FAILURE;
	render->matlevel--;
	return GZ_SUCCESS;
}

void calColor(GzRender *render, GzCoord n, GzColor &c, GzColor texture) {

	PooyanToolset::normalize(n); //normalize
	float E[] = { 0, 0, -1 };
	float s[] = { 0, 0, 0 }; //spectacular (r,g,b)
	float d[] = { 0, 0, 0 }; //diffusion (r,g,b)
	float nbackup[] = { n[0], n[1], n[2] };
	int countedLight = 0;
	for (int lCount = 0; lCount < render->numlights; lCount++) {
		// check light direction
		float nl = PooyanToolset::innerProduct(n,
				render->lights[lCount].direction);
		float ne = PooyanToolset::innerProduct(n, E);

		if (nl * ne < 0)
			continue; //Both different sign : light and eye on opposite sides of surface so that light contributes zero - skip it

		if (nl < 0 && ne < 0) { //Both negative : flip normal and compute lighting model on backside of surface
			for (int ii = 0; ii < 3; ii++)
				n[ii] *= -1;
			nl *= -1;
			ne *= -1;
		}
		//Both positive : compute lighting model
		countedLight++;
		float nl2 = 2 * nl;
		float R[] = { nl2 * n[X] - render->lights[lCount].direction[X], nl2
				* n[Y] - render->lights[lCount].direction[Y], nl2 * n[Z]
				- render->lights[lCount].direction[Z] };
		float re = PooyanToolset::innerProduct(R, E);
		if (re > 0) {
			PooyanToolset::fitValue(re, (float) 0, (float) 1);
			float t = pow(re, render->spec);
			for (int ii = 0; ii < 3; ii++)
				s[ii] += render->lights[lCount].color[ii] * t;
		}
		for (int ii = 0; ii < 3; ii++)
			d[ii] += render->lights[lCount].color[ii] * nl;
		for (int i = 0; i < 3; i++)
			n[i] = nbackup[i];
	}
	//C = (Ks sigma L [Ie (R*E)s] ) + (Kd sigma L [Ie (N*L)]) + (Ka Ia)
	for (int ii = 0; ii < 3; ii++) {
		float ka, kd, ks;
		if (render->interp_mode == GZ_NORMALS) {
			ks = render->Ks[ii];
			kd = texture[ii];
			ka = texture[ii];
		} else if (render->interp_mode == GZ_COLOR) {
			ks = 1;
			kd = 1;
			ka = 1;
		} else {
			ks = render->Ks[ii];
			kd = render->Kd[ii];
			ka = render->Ka[ii];
		}
		float kss = ks * s[ii];
		//PooyanToolset::fitValue(kss,float(0),float(1));
		float kdd = kd * d[ii];
		//PooyanToolset::fitValue(kdd,float(0),float(1));
		float kaa = ka * render->ambientlight.color[ii];
		//PooyanToolset::fitValue(kaa,float(0),float(1));
		c[ii] = kss + kdd + kaa;
		PooyanToolset::fitValue(c[ii], (float) 0, (float) 1);
	}

}

int GzDuplicateRender(GzRender *src, GzRender **dst) {
	int status = GZ_SUCCESS ;
	if (!src)
		return GZ_FAILURE;
	(*dst) = new GzRender();
	PooyanToolset::copyFromVoid(**dst, src);
	status |= GzNewDisplay(&((*dst)->display), src->display->dispClass, src->display->xres,
			src->display->yres);
	status |= GzInitDisplay((*dst)->display) ;
	//PooyanToolset::copyFromVoid((*dst)->lights,&src->lights);
	return status ;
}
