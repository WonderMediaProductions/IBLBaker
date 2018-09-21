//------------------------------------------------------------------------------------//
//                                                                                    //
//    ._____________.____   __________         __                                     //
//    |   \______   \    |  \______   \_____  |  | __ ___________                     //
//    |   ||    |  _/    |   |    |  _/\__  \ |  |/ // __ \_  __ \                    //
//    |   ||    |   \    |___|    |   \ / __ \|    <\  ___/|  | \/                    //
//    |___||______  /_______ \______  /(____  /__|_ \\___  >__|                       //
//                \/        \/      \/      \/     \/    \/                           //
//                                                                                    //
//    IBLBaker is provided under the MIT License(MIT)                                 //
//    IBLBaker uses portions of other open source software.                           //
//    Please review the LICENSE file for further details.                             //
//                                                                                    //
//    Copyright(c) 2014 Matt Davidson                                                 //
//                                                                                    //
//    Permission is hereby granted, free of charge, to any person obtaining a copy    //
//    of this software and associated documentation files(the "Software"), to deal    //
//    in the Software without restriction, including without limitation the rights    //
//    to use, copy, modify, merge, publish, distribute, sublicense, and / or sell     //
//    copies of the Software, and to permit persons to whom the Software is           //
//    furnished to do so, subject to the following conditions :                       //
//                                                                                    //
//    1. Redistributions of source code must retain the above copyright notice,       //
//    this list of conditions and the following disclaimer.                           //
//    2. Redistributions in binary form must reproduce the above copyright notice,    //
//    this list of conditions and the following disclaimer in the                     //
//    documentation and / or other materials provided with the distribution.          //
//    3. Neither the name of the copyright holder nor the names of its                //
//    contributors may be used to endorse or promote products derived                 //
//    from this software without specific prior written permission.                   //
//                                                                                    //
//    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR      //
//    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,        //
//    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.IN NO EVENT SHALL THE      //
//    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER          //
//    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,   //
//    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN       //
//    THE SOFTWARE.                                                                   //
//                                                                                    //
//------------------------------------------------------------------------------------//
float ShaderName : SHADERNAME 
<
    string ToString = "ImportanceSamplingSpecular";
> = 0.0;


float4x4 mWorldViewProj : WORLDVIEWPROJECTION;
float4x4 mWorld : WORLD;
float4x4 mView : VIEW;
float4x4 mProj : PROJECTION;
float3 vEye : EYELOCATION;

TextureCube ConvolutionSrc : CONVOLUTIONSRC;
TextureCube LastResult : LASTRESULT;
float4x4 ConvolutionViews[6] : CUBEVIEWS;
float MaxLod : IBLSOURCEMIPCOUNT;
float4 IBLCorrection : IBLCORRECTION;
float ConvolutionSamplesOffset = 0;
float ConvolutionSampleCount = 0;
float ConvolutionMaxSamples = 0;
float ConvolutionMip = 0;
uint2 ConvolutionRandom = 0;

float EnvironmentScale : IBLSOURCEENVIRONMENTSCALE;
float4 IblMaxValue : IBLMAXVALUE;

SamplerState EnvMapSampler
{
    Filter = ANISOTROPIC;
    MaxAnisotropy = 16;
    AddressU = Wrap;
    AddressV = Wrap;
};

SamplerState AccumSampler
{
    Filter = MIN_MAG_MIP_POINT;
    AddressU = Clamp;
    AddressV = Clamp;
};

struct VS_CUBEMAP_IN
{
    float4 Pos      : POSITION;
    float3 Normal   : NORMAL;
    float2 Tex      : TEXCOORD0;
};

struct GS_CUBEMAP_IN
{
    float4 Pos      : SV_POSITION;
    float3 Normal   : NORMAL;
    float2 Tex      : TEXCOORD0;
};


struct PS_CUBEMAP_IN
{
    float4 Pos      : SV_POSITION;  
    float3 Normal   : NORMAL;
    float2 Tex      : TEXCOORD0; 
    uint   Face     : TEXCOORD1;
    uint RTIndex    : SV_RenderTargetArrayIndex;
};

float3 fix_cube_lookup(float3 v, float scale) {
   float M = max(max(abs(v.x), abs(v.y)), abs(v.z));
   if (abs(v.x) != M) v.x *= scale;
   if (abs(v.y) != M) v.y *= scale;
   if (abs(v.z) != M) v.z *= scale;
   return v;
}

GS_CUBEMAP_IN VS_CubeMap( VS_CUBEMAP_IN input )
{
    GS_CUBEMAP_IN output = (GS_CUBEMAP_IN)0.0f;
    output.Pos = mul( input.Pos, mWorld );
    output.Normal = normalize(mul(input.Normal, (float3x3)(mWorld)).xyz);
    output.Tex = input.Tex;

    return output;
}

[maxvertexcount(18)]
void GS_CubeMap( triangle GS_CUBEMAP_IN input[3], inout TriangleStream<PS_CUBEMAP_IN> CubeMapStream )
{
    for (int face = 0; face < 6; face++)
    {
        PS_CUBEMAP_IN output;
        output.RTIndex = face;
        for( int v = 0; v < 3; v++ )
        {
            output.Normal = (input[v].Normal.xyz);
            output.Pos = mul(input[v].Pos, ConvolutionViews[face]);
            output.Pos = mul( output.Pos, mProj );
            output.Tex = input[v].Tex;
            output.Face = face;
            CubeMapStream.Append( output );
        }
        CubeMapStream.RestartStrip();
    }
}

float2 hammersley_seq(uint i, uint N, uint2 random)
{
    const f = 1.0 / 0x10000;
    float x = frac((float(i) + (random.x & 0xFFFF) * f));
    float y = float(reversebits(i) ^ random.y) * 2.3283064365386963e-10;
    return float2(x, y);
}

float3x3 QuaternionToMatrix(float4 quat)
{
    float3 cross = quat.yzx * quat.zxy;
    float3 square= quat.xyz * quat.xyz;
    float3 wimag = quat.w * quat.xyz;

    square = square.xyz + square.yzx;

    float3 diag = 0.5 - square;
    float3 a = (cross + wimag);
    float3 b = (cross - wimag);

    return float3x3(
    2.0 * float3(diag.x, b.z, a.y),
    2.0 * float3(a.z, diag.y, b.x),
    2.0 * float3(b.y, a.x, diag.z));
}

// float3 rescaleHDR(float3 hdrPixel)
// {
//    if (hdrPixel.x < 0)
//     hdrPixel.x = 0;
//    if (hdrPixel.y < 0)
//     hdrPixel.y = 0;
//    if (hdrPixel.z < 0)
//     hdrPixel.z = 0;

//    float intensity  = float(dot(hdrPixel, float3(0.299f,0.587f,0.114f)));

//    if (intensity > 1)
//    {
//        hdrPixel = hdrPixel - IBLCorrection.x * (hdrPixel - IblMaxValue.rgb) * hdrPixel * (hdrPixel - (IblMaxValue.rgb * 0.5));
//    }
//    else
//    {
//        hdrPixel = hdrPixel - IBLCorrection.x * (hdrPixel - IblMaxValue.rgb) * hdrPixel * (hdrPixel - 0.5);
//    }


//    // Saturation adjustment
//    hdrPixel = lerp(intensity.xxx, hdrPixel, IBLCorrection.y);

//    // Hue adjustment      
//    const float3 root = float3(0.57735, 0.57735, 0.57735);
//    float half_angle = 0.5 * radians(IBLCorrection.z); // Hue is radians of 0 tp 360 degree
//    float4 rot_quat = float4( (root * sin(half_angle)), cos(half_angle));
//    float3x3 rot_Matrix = QuaternionToMatrix(rot_quat);     
//    hdrPixel = mul(rot_Matrix, hdrPixel);
//    hdrPixel = hdrPixel * EnvironmentScale;

//    return hdrPixel; 
// }

float3 ToColor(float3 v)
{
    return float3(v*0.5+0.5);
}

float3x3 CoordinateFrameX( float3 normal )
{
	float3 e1 = float3(0,0,1);
	float3 e2 = normalize( cross( e1, normal ) );
	float3 e3 = cross( normal, e2 );
	return float3x3( e2, e3, normal );
}

float3x3 CoordinateFrameY( float3 normal )
{
	float3 e1 = float3(1,0,0);
	float3 e2 = normalize( cross( e1, normal ) );
	float3 e3 = cross( normal, e2 );
	return float3x3( e2, e3, normal );
}

float3x3 CoordinateFrameZ( float3 normal )
{
	float3 e1 = float3(0,1,0);
	float3 e2 = normalize( cross( e1, normal ) );
	float3 e3 = cross( normal, e2 );
	return float3x3( e2, e3, normal );
}

float3 importanceSampleDiffuse(float2 Xi)
{
    float CosTheta = 1.0-Xi.y;
    float SinTheta = sqrt(1.0-CosTheta*CosTheta);
    float Phi = 2*PI*Xi.x;

    float3 H;
    H.x = SinTheta * cos( Phi );
    H.y = SinTheta * sin( Phi );
    H.z = CosTheta;

    return H;
}

float3 ImportanceSample (float3 N)
{
    const float T = 0.85;

    float3 V = N;

    float3x3 cf1 = CoordinateFrameX(N);
    float3x3 cf2 = CoordinateFrameY(N);
    float3x3 cf3 = CoordinateFrameZ(N);
    float Nz = abs(N.z);
    float ts = max(0, Nz - T) / 0.1;

    float4 result = float4(0,0,0,0);
    uint sampleId = ConvolutionSamplesOffset;
    uint cubeWidth, cubeHeight;
    ConvolutionSrc.GetDimensions(cubeWidth, cubeHeight);

    uint count = ConvolutionSampleCount;
    uint2 random = ConvolutionRandom;
    float solidAngleSample = 1.0 / (count * pdf);

    for(uint i = 0; i < count; i++ )
    {
        float2 Xi = hammersley_seq(i, count, random);
        float3 L = importanceSampleDiffuse(Xi);
        float pdf = L.z / PI;
        float solidAngleTexel = 4 * PI / (6 * cubeWidth * cubeWidth);
        float lod = 0.5 * log2((float)(solidAngleSample / solidAngleTexel));

        float3 H1 = mul(L.xyz, cf1);
        float3 c1 = ConvolutionSrc.SampleLevel(EnvMapSampler, H1, lod).rgb;

        [branch]
        if (abs(N.z) > T) 
        {
            float3 H2 = mul(L.xyz, cf2);
            float3 H3 = mul(L.xyz, cf3);
            float3 c2 = ConvolutionSrc.SampleLevel(EnvMapSampler, H2, lod).rgb;
            float3 c3 = ConvolutionSrc.SampleLevel(EnvMapSampler, H3, lod).rgb;
            float3 c = lerp(c1, (c1+c2+c3)/3, ts);
            result += float4(c, 1);
        }
        else 
        {
            result += float4(c1, 1);
        }
   }

   return (result.xyz / result.w);
}

float2 cubeNormalToUV(float3 direction, uint face)
{
    switch (face)
    {
        case 0:
            return 0.5 + 0.5 * float2(-direction.z, -direction.y) / direction.x;
        case 1:
            return 0.5 + 0.5 * float2(direction.z, -direction.y) / -direction.x;
        case 2:
            return 0.5 + 0.5 * float2(direction.x, direction.z) / direction.y;
        case 3:
            return 0.5 + 0.5 * float2(direction.x, -direction.z) / -direction.y;
        case 4:
            return 0.5 + 0.5 * float2(direction.x, -direction.y) / direction.z;
        case 5:
            return 0.5 + 0.5 * float2(-direction.x, -direction.y) / -direction.z;
        default:
            return float2(0,0);
    }
}

float4 PS_CubeMap( PS_CUBEMAP_IN input) : SV_Target
{
    float3 R = normalize(input.Normal);

    // Sample source cubemap at specified mip.
    float3 importanceSampled = ImportanceSample(R);

    float4 sampledColor = 1;

    if (ConvolutionSamplesOffset > 1e-6)
    {
        float3 lastResult = LastResult.SampleLevel(EnvMapSampler, R, 0).rgb;
        sampledColor.rgb = lerp(lastResult.xyz, importanceSampled.xyz, 1.0 / (ConvolutionSamplesOffset));
    }
    else 
    {
        sampledColor.xyz = importanceSampled.xyz;
    }

    return sampledColor;
}

//-----------
// Techniques
//-----------
technique11 basic
{
    pass p0
    {
        SetVertexShader( CompileShader( vs_4_0, VS_CubeMap() ) );
        SetGeometryShader( CompileShader( gs_4_0, GS_CubeMap() ) );
        SetPixelShader( CompileShader( ps_4_0, PS_CubeMap() ) );
    }
};
