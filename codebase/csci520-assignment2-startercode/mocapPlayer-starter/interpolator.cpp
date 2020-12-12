#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include "motion.h"
#include "interpolator.h"
#include "types.h"
#include <ctime>

Interpolator::Interpolator() {
	//Set default interpolation type
	m_InterpolationType = LINEAR;

	//set default angle representation to use for interpolation
	m_AngleRepresentation = EULER;
}

Interpolator::~Interpolator() {}

//Create interpolated motion
void Interpolator::Interpolate(Motion * pInputMotion, Motion ** pOutputMotion, int N) 
{
	//Allocate new motion
	*pOutputMotion = new Motion(pInputMotion->GetNumFrames(), pInputMotion->GetSkeleton()); 

	//Perform the interpolation
	if ((m_InterpolationType == LINEAR) && (m_AngleRepresentation == EULER))
		LinearInterpolationEuler(pInputMotion, *pOutputMotion, N);
	else if ((m_InterpolationType == LINEAR) && (m_AngleRepresentation == QUATERNION))
		LinearInterpolationQuaternion(pInputMotion, *pOutputMotion, N);
	else if ((m_InterpolationType == BEZIER) && (m_AngleRepresentation == EULER))
		BezierInterpolationEuler(pInputMotion, *pOutputMotion, N);
	else if ((m_InterpolationType == BEZIER) && (m_AngleRepresentation == QUATERNION))
		BezierInterpolationQuaternion(pInputMotion, *pOutputMotion, N);
	else
	{
		printf("Error: unknown interpolation / angle representation type.\n");
		exit(1);
	}
}

// Matrix multiplication for 3 X 3 Matrix
void Interpolator::MatrixMultiplication(matrix3x3Double inputA, matrix3x3Double inputB, matrix3x3Double outputC) {
	for (int i = 0; i <= 2; i++)
		for (int j = 0; j <= 2; j++)
			outputC[i][j] = 0;

	for (int i = 0; i <= 2; i++)
		for (int j = 0; j <= 2; j++)
			for (int k = 0; k <= 2; k++)
				outputC[i][j] += (inputA[i][k] * inputB[k][j]);
}

// Interpolation type: Linear
// Angle representation for interpolation: Euler
void Interpolator::LinearInterpolationEuler(Motion * pInputMotion, Motion * pOutputMotion, int N) {
	int inputLength = pInputMotion->GetNumFrames(); // frames are indexed 0, ..., inputLength-1

	int startKeyframe = 0;

	// timer:
	clock_t startClock = clock();

	while (startKeyframe + N + 1 < inputLength)
	{
		int endKeyframe = startKeyframe + N + 1;

		Posture * startPosture = pInputMotion->GetPosture(startKeyframe);
		Posture * endPosture = pInputMotion->GetPosture(endKeyframe);

		// copy start and end keyframe
		pOutputMotion->SetPosture(startKeyframe, *startPosture);
		pOutputMotion->SetPosture(endKeyframe, *endPosture);

		// interpolate in between
		for(int frame=1; frame<=N; frame++)
		{
			Posture interpolatedPosture;
			double t = 1.0 * frame / (N+1);

			// interpolate root position
			interpolatedPosture.root_pos = startPosture->root_pos * (1-t) + endPosture->root_pos * t;

			// interpolate bone rotations
			for (int bone = 0; bone < MAX_BONES_IN_ASF_FILE; bone++)
				interpolatedPosture.bone_rotation[bone] = startPosture->bone_rotation[bone] * (1-t) + endPosture->bone_rotation[bone] * t;

			pOutputMotion->SetPosture(startKeyframe + frame, interpolatedPosture);
		}

		startKeyframe = endKeyframe;
	}

	for(int frame=startKeyframe+1; frame<inputLength; frame++)
		pOutputMotion->SetPosture(frame, *(pInputMotion->GetPosture(frame)));

	double timeTaken = clock() - startClock;
	printf("Computation Time: %lf\n", timeTaken / (double) CLOCKS_PER_SEC);
}

// Interpolation type: Linear
// Angle representation for interpolation: Quaternion
void Interpolator::LinearInterpolationQuaternion(Motion * pInputMotion, Motion * pOutputMotion, int N) {
	int numOfFrames;
	int firstFrame, lastFrame;
	double t;
	Quaternion<double> firstFrameQuaternion, lastFrameQuaternion;

	// Get total number of frames
	numOfFrames = pInputMotion->GetNumFrames();

	// Start the timer
	clock_t startClock = clock();

	for (firstFrame = 0; firstFrame + N + 1 < numOfFrames; firstFrame = lastFrame){
		lastFrame = firstFrame + N + 1;

		Posture * firstFramePosture = pInputMotion->GetPosture(firstFrame);
		pOutputMotion->SetPosture(firstFrame, *firstFramePosture);

		Posture * lastFramePosture = pInputMotion->GetPosture(lastFrame);
		pOutputMotion->SetPosture(lastFrame, *lastFramePosture);

		// Frame interpolation loop
		for (int frameIndex = 1; frameIndex <= N; frameIndex++) {
			Posture posture;

			// Interpolate root position
			t = 1.0 * frameIndex / (N + 1);
			posture.root_pos = firstFramePosture->root_pos * (1 - t) + lastFramePosture->root_pos * t;

			// Bone interpolation loop
			for (int boneIndex = 0; boneIndex < MAX_BONES_IN_ASF_FILE; boneIndex++) { 
				// Convert Euler to Quaternion
				Euler2Quaternion(firstFramePosture->bone_rotation[boneIndex].p, firstFrameQuaternion);
				Euler2Quaternion(lastFramePosture->bone_rotation[boneIndex].p, lastFrameQuaternion);

				// Slerp
				Quaternion<double> slerp;
				slerp = Slerp(t, firstFrameQuaternion, lastFrameQuaternion);

				// Convert Quaternion to Euler and set new rotation values
				Quaternion2Euler(slerp, posture.bone_rotation[boneIndex].p);
			}
			pOutputMotion->SetPosture(firstFrame + frameIndex, posture);
		}
	}

	for (int frameIndex = firstFrame + 1; frameIndex < numOfFrames; frameIndex++)
		pOutputMotion->SetPosture(frameIndex, *(pInputMotion->GetPosture(frameIndex)));

	double timeTaken = clock() - startClock;
	printf("Computation Time: %lf\n", timeTaken / (double) CLOCKS_PER_SEC);
}

// Interpolation type: Bezier
// Angle representation for interpolation: Euler
void Interpolator::BezierInterpolationEuler(Motion * pInputMotion, Motion * pOutputMotion, int N) {
	int numOfFrames;
	int firstFrame, lastFrame, prevFrame, nextFrame;

	double t;

	Quaternion<double> firstFrameQuaternion, lastFrameQuaternion;
	
	vector firstFrameEuler, lastFrameEuler, prevFrameEuler, nextFrameEuler, tempEuler;

	// Get total number of frames
	numOfFrames = pInputMotion->GetNumFrames();

	// Start the timer
	clock_t startClock = clock();

	for(firstFrame = 0; firstFrame + N + 1 < numOfFrames; firstFrame = lastFrame) {
		lastFrame = firstFrame + N + 1;

		prevFrame = firstFrame - N - 1;
		if (firstFrame == 0) {
			prevFrame = firstFrame;
		}

		vector * previousFrameRotation = pInputMotion->GetPosture(prevFrame)->bone_rotation;

		nextFrame = lastFrame + N + 1;
		if (lastFrame + N + 1 >= numOfFrames) {
			nextFrame = lastFrame;
		}

		vector * nextFrameRotation = pInputMotion->GetPosture(nextFrame)->bone_rotation;

		Posture * firstFramePosture = pInputMotion->GetPosture(firstFrame);
		pOutputMotion->SetPosture(firstFrame, *firstFramePosture);

		Posture * lastFramePosture = pInputMotion->GetPosture(lastFrame);
		pOutputMotion->SetPosture(lastFrame, *lastFramePosture);

		// Frame interpolation loop
		for (int frameIndex = 1; frameIndex <= N; frameIndex++) {
			Posture posture;

			// Interpolate root position
			t = 1.0 * frameIndex / (N + 1);
			posture.root_pos = firstFramePosture->root_pos * (1 - t) + lastFramePosture->root_pos * t;

			// Bone interpolation loop
			for (int boneIndex = 0; boneIndex < MAX_BONES_IN_ASF_FILE; boneIndex++) {
				prevFrameEuler = previousFrameRotation[boneIndex];
				nextFrameEuler = nextFrameRotation[boneIndex];

				firstFrameEuler = firstFramePosture->bone_rotation[boneIndex];
				lastFrameEuler = lastFramePosture->bone_rotation[boneIndex];

				prevFrameEuler = firstFrameEuler + firstFrameEuler - prevFrameEuler;
				prevFrameEuler = (prevFrameEuler + lastFrameEuler) * 0.5;
				prevFrameEuler = firstFrameEuler + (prevFrameEuler - firstFrameEuler) / 3.0;
				
				tempEuler = lastFrameEuler + lastFrameEuler - firstFrameEuler;
				nextFrameEuler = (tempEuler + nextFrameEuler) * 0.5;
				nextFrameEuler = lastFrameEuler + lastFrameEuler - nextFrameEuler;
				nextFrameEuler = lastFrameEuler + (nextFrameEuler - lastFrameEuler) / 3.0;

				posture.bone_rotation[boneIndex] = DeCasteljauEuler(t, firstFrameEuler, prevFrameEuler, nextFrameEuler, lastFrameEuler);
			}
			pOutputMotion->SetPosture(firstFrame + frameIndex, posture);
		}
	}
	
	for (int frameIndex = firstFrame + 1; frameIndex < numOfFrames; frameIndex++)
		pOutputMotion->SetPosture(frameIndex, *(pInputMotion->GetPosture(frameIndex)));

	double timeTaken = clock() - startClock;
	printf("Computation Time: %lf\n", timeTaken / (double) CLOCKS_PER_SEC);
}

// Interpolation type: Bezier
// Angle representation for interpolation: Quaternion
void Interpolator::BezierInterpolationQuaternion(Motion * pInputMotion, Motion * pOutputMotion, int N) {
	int numOfFrames;
	int firstFrame, lastFrame, previousFrame, nextFrame;

	double t;
	double firstFrameAngles[3], lastFrameAngles[3], previousFrameAngles[3], nextFrameAngles[3];
	double finalEulerAngle[3];

	Quaternion<double> firstFrameQuaternion, lastFrameQuaternion, previousFrameQuaternion, nextFrameQuaternion, tempQuaternion;
	Quaternion<double> finalQuaternionAngle;

	vector first, last;

	// Get total number of frames
	numOfFrames = pInputMotion->GetNumFrames();

	// Start the timer
	clock_t startClock = clock();

	for (firstFrame = 0; firstFrame + N + 1 < numOfFrames; firstFrame = lastFrame) {
		lastFrame = firstFrame + N + 1;

		previousFrame = firstFrame - N - 1;
		if (firstFrame == 0) {
			previousFrame = firstFrame;
		}

		vector *previousFrameRotation = pInputMotion->GetPosture(previousFrame)->bone_rotation;

		nextFrame = lastFrame + N + 1;
		if (lastFrame + N + 1 >= numOfFrames) {
			nextFrame = lastFrame;
		}

		vector *nextFrameRotation = pInputMotion->GetPosture(nextFrame)->bone_rotation;

		Posture * firstFramePosture = pInputMotion->GetPosture(firstFrame);
		pOutputMotion->SetPosture(firstFrame, *firstFramePosture);

		Posture * lastFramePosture = pInputMotion->GetPosture(lastFrame);
		pOutputMotion->SetPosture(lastFrame, *lastFramePosture);

		// Frame interpolation loop
		for (int frameIndex = 1; frameIndex <= N; frameIndex++){
			Posture posture;

			// Interpolate root position
			t = 1.0 * frameIndex / (N + 1);
			posture.root_pos = firstFramePosture->root_pos * (1 - t) + lastFramePosture->root_pos * t;

			// Bone interpolation loop
			for (int boneIndex = 0; boneIndex < MAX_BONES_IN_ASF_FILE; boneIndex++){
				// Get first, last, previous and next frame Euler angles and convert them to Quaternions
				firstFramePosture->bone_rotation[boneIndex].getValue(firstFrameAngles);
				Euler2Quaternion(firstFrameAngles, firstFrameQuaternion);

				lastFramePosture->bone_rotation[boneIndex].getValue(lastFrameAngles);
				Euler2Quaternion(lastFrameAngles, lastFrameQuaternion);

				previousFrameRotation[boneIndex].getValue(previousFrameAngles);
				Euler2Quaternion(previousFrameAngles, previousFrameQuaternion);

				nextFrameRotation[boneIndex].getValue(nextFrameAngles);
				Euler2Quaternion(nextFrameAngles, nextFrameQuaternion);

				previousFrameQuaternion = Double(previousFrameQuaternion, firstFrameQuaternion);
				previousFrameQuaternion = Slerp(0.5, previousFrameQuaternion, lastFrameQuaternion);
				previousFrameQuaternion = Slerp(1.0 / 3.0, firstFrameQuaternion, previousFrameQuaternion);

				tempQuaternion = Double(firstFrameQuaternion, lastFrameQuaternion);
				nextFrameQuaternion = Slerp(0.5, tempQuaternion, nextFrameQuaternion);
				nextFrameQuaternion = Double(nextFrameQuaternion, lastFrameQuaternion);
				nextFrameQuaternion = Slerp(1.0 / 3.0, lastFrameQuaternion, nextFrameQuaternion);

				finalQuaternionAngle = DeCasteljauQuaternion(t, firstFrameQuaternion, previousFrameQuaternion, nextFrameQuaternion, lastFrameQuaternion);
				
				// Convert Quaternion to Euler and set new rotation values
				Quaternion2Euler(finalQuaternionAngle, finalEulerAngle);
				posture.bone_rotation[boneIndex].setValue(finalEulerAngle);
			}
			pOutputMotion->SetPosture(firstFrame + frameIndex, posture);
		}
	}

	for (int frameIndex = firstFrame + 1; frameIndex < numOfFrames; frameIndex++) {
		pOutputMotion->SetPosture(frameIndex, *(pInputMotion->GetPosture(frameIndex)));
	}

	double timeTaken = clock() - startClock;
	printf("Computation Time: %lf\n", timeTaken / (double) CLOCKS_PER_SEC);
}

// Convert Rotation Matrix to Euler
void Interpolator::Rotation2Euler(double R[9], double angles[3]) {
	double cy = sqrt(R[0]*R[0] + R[3]*R[3]);

	if (cy > 16*DBL_EPSILON) {
		angles[0] = atan2(R[7], R[8]);
		angles[1] = atan2(-R[6], cy);
		angles[2] = atan2(R[3], R[0]);
	} 
	else {
		angles[0] = atan2(-R[5], R[4]);
		angles[1] = atan2(-R[6], cy);
		angles[2] = 0;
	}

	for(int i=0; i<3; i++)
		angles[i] *= 180 / M_PI;
}

// Convert Euler to Rotation Matrix
void Interpolator::Euler2Rotation(double angles[3], double R[9]) {
	int i, j;

	double radians[3] = { ConvertDegreeToRadian(angles[0]), ConvertDegreeToRadian(angles[1]), ConvertDegreeToRadian(angles[2]) };
	matrix3x3Double tempOutput, finalOutput;

	matrix3x3Double M1 = { 
		{ cos(radians[2]), 	-sin(radians[2]), 	0.0 },
		{ sin(radians[2]), 	cos(radians[2]), 	0.0 }, 
		{ 0.0, 				0.0, 				1.0 }
	};
	matrix3x3Double M2 = { 
		{ cos(radians[1]), 	0.0, 	sin(radians[1]) }, 
		{ 0.0, 				1.0, 	0.0 }, 
		{-sin(radians[1]), 	0.0, 	cos(radians[1]) }
	};

	MatrixMultiplication(M1, M2, tempOutput);

	matrix3x3Double M3 = { 
		{	1.0, 	0.0, 				0.0 }, 
		{ 	0.0, 	cos(radians[0]), 	-sin(radians[0]) }, 
		{	0.0, 	sin(radians[0]), 	cos(radians[0]) }
	};

	MatrixMultiplication(tempOutput, M3, finalOutput);

	for (i = 0; i <= 2; i++)
		for (j = 0; j <= 2; j++)
			R[j + (3*i)] = finalOutput[i][j];
}

// Covert Euler to Quaternion
void Interpolator::Euler2Quaternion(double angles[3], Quaternion<double> & q) {
	// Covert Euler to Rotation
	double t_matrix[9];
	Euler2Rotation(angles, t_matrix);

	// Then convert it to Quaternion
	q = Quaternion<double>::Matrix2Quaternion(t_matrix);
	q.Normalize();
}

// Covert Quaternion to Euler
void Interpolator::Quaternion2Euler(Quaternion<double> & q, double angles[3]) {
	// Quaternion to Rotation
	double t_matrix[9];
	q.Quaternion2Matrix(t_matrix);

	// Then convert to Euler.
	Rotation2Euler(t_matrix, angles);
}

// Calculate Slerp
Quaternion<double> Interpolator::Slerp(double t, Quaternion<double> & qStart, Quaternion<double> & qEnd_){
	Quaternion<double> result, newQEnd_;

	// cos(theta) = q1 . q2
	double dotProduct;
	dotProduct = qStart.Gets() * qEnd_.Gets() + qStart.Getx() * qEnd_.Getx() + qStart.Gety() * qEnd_.Gety() + qStart.Getz() * qEnd_.Getz();

	// Recalculate newQEnd_ depending on whether dotProduct is >= 0 or not
	if(dotProduct >= 0) {
		newQEnd_ = qEnd_;
	}
	else {
		dotProduct *= -1;
		newQEnd_.Set(-qEnd_.Gets(), -qEnd_.Getx(), -qEnd_.Gety(), -qEnd_.Getz());
	}

	// Find theta
	float theta;
	theta = acosf(dotProduct);

	// Calculate slerp
	if(theta != 0) {
		result = ((sinf((1 - t) * theta) / sinf(theta)) * qStart + (sinf(t * theta) / sinf(theta)) * newQEnd_);
		result.Normalize();
		return result;
	}
	else {
		return qEnd_;
	}
}

// Calculate 2(p . q)q – p where p . q is the dot product
Quaternion<double> Interpolator::Double(Quaternion<double> p, Quaternion<double> q) {
	double dotProduct = p.Gets() * q.Gets() + p.Getx() * q.Getx() + p.Gety() * q.Gety() + p.Getz() * q.Getz();

	Quaternion<double> result;
	result = 2 * (dotProduct) * q - p;

	return result;
}

// Calculate DeCasteljau Euler : 
//     Q0 = P1 * t + P0 * (1 - t);
//     Q1 = P2 * t + P1 * (1 - t);
//     Q2 = P3 * t + P2 * (1 - t);
//     R0 = Q1 * t + Q0 * (1 - t);
//     R1 = Q2 * t + Q1 * (1 - t);
// Return: 
//     P(t) = R1 * t + R0 * (1 - t);
vector Interpolator::DeCasteljauEuler(double t, vector p0, vector p1, vector p2, vector p3) {
	vector q0 = p1 * t + p0 * (1 - t);
	vector q1 = p2 * t + p1 * (1 - t);
	vector q2 = p3 * t + p2 * (1 - t);

	vector r0 = q1 * t + q0 * (1 - t);
	vector r1 = q2 * t + q1 * (1 - t);
	
	vector result = r1 * t + r0 * (1 - t);
	return result;
}

// Calculate DeCasteljau Quaternion : 
//     Q0 = Slerp(P0,P1,t)
//     Q1 = Slerp(P1,P2,t) 
//     Q2 = Slerp(P2,P3,t) 
//     R0 = Slerp(Q0,Q1,t)
//     R1 = Slerp(Q1,Q2,t)
// Return: 
//     P(t) = Slerp(R0,R1,t)
Quaternion<double> Interpolator::DeCasteljauQuaternion(double t, Quaternion<double> p0, Quaternion<double> p1, Quaternion<double> p2, Quaternion<double> p3) {
	Quaternion<double> q0 = Slerp(t, p0, p1);
	Quaternion<double> q1 = Slerp(t, p1, p2);
	Quaternion<double> q2 = Slerp(t, p2, p3);

	Quaternion<double> r0 = Slerp(t, q0, q1);
	Quaternion<double> r1 = Slerp(t, q1, q2);

	Quaternion<double> result = Slerp(t, r0, r1);
	return result;
}

