// Copyright (c) 2015 Alberto Otero de la Roza
// <aoterodelaroza@gmail.com>,
// Ángel Martín Pendás <angel@fluor.quimica.uniovi.es> and Víctor Luaña
// <victor@fluor.quimica.uniovi.es>.
//
// critic2 is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at
// your option) any later version.
//
// critic2 is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see
// <http://www.gnu.org/licenses/>.

// Preprocessor-constants for Fortran.

#ifndef CONSTANTS_H
#define CONSTANTS_H

// threads
/* Function return values */
extern "C" const int const_thrd_error;
extern "C" const int const_thrd_success;
extern "C" const int const_thrd_timedout;
extern "C" const int const_thrd_busy;
extern "C" const int const_thrd_nomem;
/* Mutex types */
extern "C" const int const_mtx_plain;
extern "C" const int const_mtx_timed;
extern "C" const int const_mtx_recursive;

// OpenGL
extern "C" const int const_GL_DEPTH_BUFFER_BIT;
extern "C" const int const_GL_STENCIL_BUFFER_BIT;
extern "C" const int const_GL_COLOR_BUFFER_BIT;
extern "C" const int const_GL_FALSE;
extern "C" const int const_GL_TRUE;
extern "C" const int const_GL_POINTS;
extern "C" const int const_GL_LINES;
extern "C" const int const_GL_LINE_LOOP;
extern "C" const int const_GL_LINE_STRIP;
extern "C" const int const_GL_TRIANGLES;
extern "C" const int const_GL_TRIANGLE_STRIP;
extern "C" const int const_GL_TRIANGLE_FAN;
extern "C" const int const_GL_QUADS;
extern "C" const int const_GL_NEVER;
extern "C" const int const_GL_LESS;
extern "C" const int const_GL_EQUAL;
extern "C" const int const_GL_LEQUAL;
extern "C" const int const_GL_GREATER;
extern "C" const int const_GL_NOTEQUAL;
extern "C" const int const_GL_GEQUAL;
extern "C" const int const_GL_ALWAYS;
extern "C" const int const_GL_ZERO;
extern "C" const int const_GL_ONE;
extern "C" const int const_GL_SRC_COLOR;
extern "C" const int const_GL_ONE_MINUS_SRC_COLOR;
extern "C" const int const_GL_SRC_ALPHA;
extern "C" const int const_GL_ONE_MINUS_SRC_ALPHA;
extern "C" const int const_GL_DST_ALPHA;
extern "C" const int const_GL_ONE_MINUS_DST_ALPHA;
extern "C" const int const_GL_DST_COLOR;
extern "C" const int const_GL_ONE_MINUS_DST_COLOR;
extern "C" const int const_GL_SRC_ALPHA_SATURATE;
extern "C" const int const_GL_NONE;
extern "C" const int const_GL_FRONT_LEFT;
extern "C" const int const_GL_FRONT_RIGHT;
extern "C" const int const_GL_BACK_LEFT;
extern "C" const int const_GL_BACK_RIGHT;
extern "C" const int const_GL_FRONT;
extern "C" const int const_GL_BACK;
extern "C" const int const_GL_LEFT;
extern "C" const int const_GL_RIGHT;
extern "C" const int const_GL_FRONT_AND_BACK;
extern "C" const int const_GL_NO_ERROR;
extern "C" const int const_GL_INVALID_ENUM;
extern "C" const int const_GL_INVALID_VALUE;
extern "C" const int const_GL_INVALID_OPERATION;
extern "C" const int const_GL_OUT_OF_MEMORY;
extern "C" const int const_GL_CW;
extern "C" const int const_GL_CCW;
extern "C" const int const_GL_POINT_SIZE;
extern "C" const int const_GL_POINT_SIZE_RANGE;
extern "C" const int const_GL_POINT_SIZE_GRANULARITY;
extern "C" const int const_GL_LINE_SMOOTH;
extern "C" const int const_GL_LINE_WIDTH;
extern "C" const int const_GL_LINE_WIDTH_RANGE;
extern "C" const int const_GL_LINE_WIDTH_GRANULARITY;
extern "C" const int const_GL_POLYGON_MODE;
extern "C" const int const_GL_POLYGON_SMOOTH;
extern "C" const int const_GL_CULL_FACE;
extern "C" const int const_GL_CULL_FACE_MODE;
extern "C" const int const_GL_FRONT_FACE;
extern "C" const int const_GL_DEPTH_RANGE;
extern "C" const int const_GL_DEPTH_TEST;
extern "C" const int const_GL_DEPTH_WRITEMASK;
extern "C" const int const_GL_DEPTH_CLEAR_VALUE;
extern "C" const int const_GL_DEPTH_FUNC;
extern "C" const int const_GL_STENCIL_TEST;
extern "C" const int const_GL_STENCIL_CLEAR_VALUE;
extern "C" const int const_GL_STENCIL_FUNC;
extern "C" const int const_GL_STENCIL_VALUE_MASK;
extern "C" const int const_GL_STENCIL_FAIL;
extern "C" const int const_GL_STENCIL_PASS_DEPTH_FAIL;
extern "C" const int const_GL_STENCIL_PASS_DEPTH_PASS;
extern "C" const int const_GL_STENCIL_REF;
extern "C" const int const_GL_STENCIL_WRITEMASK;
extern "C" const int const_GL_VIEWPORT;
extern "C" const int const_GL_DITHER;
extern "C" const int const_GL_BLEND_DST;
extern "C" const int const_GL_BLEND_SRC;
extern "C" const int const_GL_BLEND;
extern "C" const int const_GL_LOGIC_OP_MODE;
extern "C" const int const_GL_COLOR_LOGIC_OP;
extern "C" const int const_GL_DRAW_BUFFER;
extern "C" const int const_GL_READ_BUFFER;
extern "C" const int const_GL_SCISSOR_BOX;
extern "C" const int const_GL_SCISSOR_TEST;
extern "C" const int const_GL_COLOR_CLEAR_VALUE;
extern "C" const int const_GL_COLOR_WRITEMASK;
extern "C" const int const_GL_DOUBLEBUFFER;
extern "C" const int const_GL_STEREO;
extern "C" const int const_GL_LINE_SMOOTH_HINT;
extern "C" const int const_GL_POLYGON_SMOOTH_HINT;
extern "C" const int const_GL_UNPACK_SWAP_BYTES;
extern "C" const int const_GL_UNPACK_LSB_FIRST;
extern "C" const int const_GL_UNPACK_ROW_LENGTH;
extern "C" const int const_GL_UNPACK_SKIP_ROWS;
extern "C" const int const_GL_UNPACK_SKIP_PIXELS;
extern "C" const int const_GL_UNPACK_ALIGNMENT;
extern "C" const int const_GL_PACK_SWAP_BYTES;
extern "C" const int const_GL_PACK_LSB_FIRST;
extern "C" const int const_GL_PACK_ROW_LENGTH;
extern "C" const int const_GL_PACK_SKIP_ROWS;
extern "C" const int const_GL_PACK_SKIP_PIXELS;
extern "C" const int const_GL_PACK_ALIGNMENT;
extern "C" const int const_GL_MAX_TEXTURE_SIZE;
extern "C" const int const_GL_MAX_VIEWPORT_DIMS;
extern "C" const int const_GL_SUBPIXEL_BITS;
extern "C" const int const_GL_TEXTURE_1D;
extern "C" const int const_GL_TEXTURE_2D;
extern "C" const int const_GL_POLYGON_OFFSET_UNITS;
extern "C" const int const_GL_POLYGON_OFFSET_POINT;
extern "C" const int const_GL_POLYGON_OFFSET_LINE;
extern "C" const int const_GL_POLYGON_OFFSET_FILL;
extern "C" const int const_GL_POLYGON_OFFSET_FACTOR;
extern "C" const int const_GL_TEXTURE_BINDING_1D;
extern "C" const int const_GL_TEXTURE_BINDING_2D;
extern "C" const int const_GL_TEXTURE_WIDTH;
extern "C" const int const_GL_TEXTURE_HEIGHT;
extern "C" const int const_GL_TEXTURE_INTERNAL_FORMAT;
extern "C" const int const_GL_TEXTURE_BORDER_COLOR;
extern "C" const int const_GL_TEXTURE_RED_SIZE;
extern "C" const int const_GL_TEXTURE_GREEN_SIZE;
extern "C" const int const_GL_TEXTURE_BLUE_SIZE;
extern "C" const int const_GL_TEXTURE_ALPHA_SIZE;
extern "C" const int const_GL_DONT_CARE;
extern "C" const int const_GL_FASTEST;
extern "C" const int const_GL_NICEST;
extern "C" const int const_GL_BYTE;
extern "C" const int const_GL_UNSIGNED_BYTE;
extern "C" const int const_GL_SHORT;
extern "C" const int const_GL_UNSIGNED_SHORT;
extern "C" const int const_GL_INT;
extern "C" const int const_GL_UNSIGNED_INT;
extern "C" const int const_GL_FLOAT;
extern "C" const int const_GL_DOUBLE;
extern "C" const int const_GL_STACK_OVERFLOW;
extern "C" const int const_GL_STACK_UNDERFLOW;
extern "C" const int const_GL_CLEAR;
extern "C" const int const_GL_AND;
extern "C" const int const_GL_AND_REVERSE;
extern "C" const int const_GL_COPY;
extern "C" const int const_GL_AND_INVERTED;
extern "C" const int const_GL_NOOP;
extern "C" const int const_GL_XOR;
extern "C" const int const_GL_OR;
extern "C" const int const_GL_NOR;
extern "C" const int const_GL_EQUIV;
extern "C" const int const_GL_INVERT;
extern "C" const int const_GL_OR_REVERSE;
extern "C" const int const_GL_COPY_INVERTED;
extern "C" const int const_GL_OR_INVERTED;
extern "C" const int const_GL_NAND;
extern "C" const int const_GL_SET;
extern "C" const int const_GL_TEXTURE;
extern "C" const int const_GL_COLOR;
extern "C" const int const_GL_DEPTH;
extern "C" const int const_GL_STENCIL;
extern "C" const int const_GL_STENCIL_INDEX;
extern "C" const int const_GL_DEPTH_COMPONENT;
extern "C" const int const_GL_RED;
extern "C" const int const_GL_GREEN;
extern "C" const int const_GL_BLUE;
extern "C" const int const_GL_ALPHA;
extern "C" const int const_GL_RGB;
extern "C" const int const_GL_RGBA;
extern "C" const int const_GL_POINT;
extern "C" const int const_GL_LINE;
extern "C" const int const_GL_FILL;
extern "C" const int const_GL_KEEP;
extern "C" const int const_GL_REPLACE;
extern "C" const int const_GL_INCR;
extern "C" const int const_GL_DECR;
extern "C" const int const_GL_VENDOR;
extern "C" const int const_GL_RENDERER;
extern "C" const int const_GL_VERSION;
extern "C" const int const_GL_EXTENSIONS;
extern "C" const int const_GL_NEAREST;
extern "C" const int const_GL_LINEAR;
extern "C" const int const_GL_NEAREST_MIPMAP_NEAREST;
extern "C" const int const_GL_LINEAR_MIPMAP_NEAREST;
extern "C" const int const_GL_NEAREST_MIPMAP_LINEAR;
extern "C" const int const_GL_LINEAR_MIPMAP_LINEAR;
extern "C" const int const_GL_TEXTURE_MAG_FILTER;
extern "C" const int const_GL_TEXTURE_MIN_FILTER;
extern "C" const int const_GL_TEXTURE_WRAP_S;
extern "C" const int const_GL_TEXTURE_WRAP_T;
extern "C" const int const_GL_PROXY_TEXTURE_1D;
extern "C" const int const_GL_PROXY_TEXTURE_2D;
extern "C" const int const_GL_REPEAT;
extern "C" const int const_GL_R3_G3_B2;
extern "C" const int const_GL_RGB4;
extern "C" const int const_GL_RGB5;
extern "C" const int const_GL_RGB8;
extern "C" const int const_GL_RGB10;
extern "C" const int const_GL_RGB12;
extern "C" const int const_GL_RGB16;
extern "C" const int const_GL_RGBA2;
extern "C" const int const_GL_RGBA4;
extern "C" const int const_GL_RGB5_A1;
extern "C" const int const_GL_RGBA8;
extern "C" const int const_GL_RGB10_A2;
extern "C" const int const_GL_RGBA12;
extern "C" const int const_GL_RGBA16;
extern "C" const int const_GL_UNSIGNED_BYTE_3_3_2;
extern "C" const int const_GL_UNSIGNED_SHORT_4_4_4_4;
extern "C" const int const_GL_UNSIGNED_SHORT_5_5_5_1;
extern "C" const int const_GL_UNSIGNED_INT_8_8_8_8;
extern "C" const int const_GL_UNSIGNED_INT_10_10_10_2;
extern "C" const int const_GL_TEXTURE_BINDING_3D;
extern "C" const int const_GL_PACK_SKIP_IMAGES;
extern "C" const int const_GL_PACK_IMAGE_HEIGHT;
extern "C" const int const_GL_UNPACK_SKIP_IMAGES;
extern "C" const int const_GL_UNPACK_IMAGE_HEIGHT;
extern "C" const int const_GL_TEXTURE_3D;
extern "C" const int const_GL_PROXY_TEXTURE_3D;
extern "C" const int const_GL_TEXTURE_DEPTH;
extern "C" const int const_GL_TEXTURE_WRAP_R;
extern "C" const int const_GL_MAX_3D_TEXTURE_SIZE;
extern "C" const int const_GL_UNSIGNED_BYTE_2_3_3_REV;
extern "C" const int const_GL_UNSIGNED_SHORT_5_6_5;
extern "C" const int const_GL_UNSIGNED_SHORT_5_6_5_REV;
extern "C" const int const_GL_UNSIGNED_SHORT_4_4_4_4_REV;
extern "C" const int const_GL_UNSIGNED_SHORT_1_5_5_5_REV;
extern "C" const int const_GL_UNSIGNED_INT_8_8_8_8_REV;
extern "C" const int const_GL_UNSIGNED_INT_2_10_10_10_REV;
extern "C" const int const_GL_BGR;
extern "C" const int const_GL_BGRA;
extern "C" const int const_GL_MAX_ELEMENTS_VERTICES;
extern "C" const int const_GL_MAX_ELEMENTS_INDICES;
extern "C" const int const_GL_CLAMP_TO_EDGE;
extern "C" const int const_GL_TEXTURE_MIN_LOD;
extern "C" const int const_GL_TEXTURE_MAX_LOD;
extern "C" const int const_GL_TEXTURE_BASE_LEVEL;
extern "C" const int const_GL_TEXTURE_MAX_LEVEL;
extern "C" const int const_GL_SMOOTH_POINT_SIZE_RANGE;
extern "C" const int const_GL_SMOOTH_POINT_SIZE_GRANULARITY;
extern "C" const int const_GL_SMOOTH_LINE_WIDTH_RANGE;
extern "C" const int const_GL_SMOOTH_LINE_WIDTH_GRANULARITY;
extern "C" const int const_GL_ALIASED_LINE_WIDTH_RANGE;
extern "C" const int const_GL_CONSTANT_COLOR;
extern "C" const int const_GL_ONE_MINUS_CONSTANT_COLOR;
extern "C" const int const_GL_CONSTANT_ALPHA;
extern "C" const int const_GL_ONE_MINUS_CONSTANT_ALPHA;
extern "C" const int const_GL_BLEND_COLOR;
extern "C" const int const_GL_FUNC_ADD;
extern "C" const int const_GL_MIN;
extern "C" const int const_GL_MAX;
extern "C" const int const_GL_BLEND_EQUATION;
extern "C" const int const_GL_FUNC_SUBTRACT;
extern "C" const int const_GL_FUNC_REVERSE_SUBTRACT;
extern "C" const int const_GL_TEXTURE0;
extern "C" const int const_GL_TEXTURE1;
extern "C" const int const_GL_TEXTURE2;
extern "C" const int const_GL_TEXTURE3;
extern "C" const int const_GL_TEXTURE4;
extern "C" const int const_GL_TEXTURE5;
extern "C" const int const_GL_TEXTURE6;
extern "C" const int const_GL_TEXTURE7;
extern "C" const int const_GL_TEXTURE8;
extern "C" const int const_GL_TEXTURE9;
extern "C" const int const_GL_TEXTURE10;
extern "C" const int const_GL_TEXTURE11;
extern "C" const int const_GL_TEXTURE12;
extern "C" const int const_GL_TEXTURE13;
extern "C" const int const_GL_TEXTURE14;
extern "C" const int const_GL_TEXTURE15;
extern "C" const int const_GL_TEXTURE16;
extern "C" const int const_GL_TEXTURE17;
extern "C" const int const_GL_TEXTURE18;
extern "C" const int const_GL_TEXTURE19;
extern "C" const int const_GL_TEXTURE20;
extern "C" const int const_GL_TEXTURE21;
extern "C" const int const_GL_TEXTURE22;
extern "C" const int const_GL_TEXTURE23;
extern "C" const int const_GL_TEXTURE24;
extern "C" const int const_GL_TEXTURE25;
extern "C" const int const_GL_TEXTURE26;
extern "C" const int const_GL_TEXTURE27;
extern "C" const int const_GL_TEXTURE28;
extern "C" const int const_GL_TEXTURE29;
extern "C" const int const_GL_TEXTURE30;
extern "C" const int const_GL_TEXTURE31;
extern "C" const int const_GL_ACTIVE_TEXTURE;
extern "C" const int const_GL_MULTISAMPLE;
extern "C" const int const_GL_SAMPLE_ALPHA_TO_COVERAGE;
extern "C" const int const_GL_SAMPLE_ALPHA_TO_ONE;
extern "C" const int const_GL_SAMPLE_COVERAGE;
extern "C" const int const_GL_SAMPLE_BUFFERS;
extern "C" const int const_GL_SAMPLES;
extern "C" const int const_GL_SAMPLE_COVERAGE_VALUE;
extern "C" const int const_GL_SAMPLE_COVERAGE_INVERT;
extern "C" const int const_GL_TEXTURE_CUBE_MAP;
extern "C" const int const_GL_TEXTURE_BINDING_CUBE_MAP;
extern "C" const int const_GL_TEXTURE_CUBE_MAP_POSITIVE_X;
extern "C" const int const_GL_TEXTURE_CUBE_MAP_NEGATIVE_X;
extern "C" const int const_GL_TEXTURE_CUBE_MAP_POSITIVE_Y;
extern "C" const int const_GL_TEXTURE_CUBE_MAP_NEGATIVE_Y;
extern "C" const int const_GL_TEXTURE_CUBE_MAP_POSITIVE_Z;
extern "C" const int const_GL_TEXTURE_CUBE_MAP_NEGATIVE_Z;
extern "C" const int const_GL_PROXY_TEXTURE_CUBE_MAP;
extern "C" const int const_GL_MAX_CUBE_MAP_TEXTURE_SIZE;
extern "C" const int const_GL_COMPRESSED_RGB;
extern "C" const int const_GL_COMPRESSED_RGBA;
extern "C" const int const_GL_TEXTURE_COMPRESSION_HINT;
extern "C" const int const_GL_TEXTURE_COMPRESSED_IMAGE_SIZE;
extern "C" const int const_GL_TEXTURE_COMPRESSED;
extern "C" const int const_GL_NUM_COMPRESSED_TEXTURE_FORMATS;
extern "C" const int const_GL_COMPRESSED_TEXTURE_FORMATS;
extern "C" const int const_GL_CLAMP_TO_BORDER;
extern "C" const int const_GL_BLEND_DST_RGB;
extern "C" const int const_GL_BLEND_SRC_RGB;
extern "C" const int const_GL_BLEND_DST_ALPHA;
extern "C" const int const_GL_BLEND_SRC_ALPHA;
extern "C" const int const_GL_POINT_FADE_THRESHOLD_SIZE;
extern "C" const int const_GL_DEPTH_COMPONENT16;
extern "C" const int const_GL_DEPTH_COMPONENT24;
extern "C" const int const_GL_DEPTH_COMPONENT32;
extern "C" const int const_GL_MIRRORED_REPEAT;
extern "C" const int const_GL_MAX_TEXTURE_LOD_BIAS;
extern "C" const int const_GL_TEXTURE_LOD_BIAS;
extern "C" const int const_GL_INCR_WRAP;
extern "C" const int const_GL_DECR_WRAP;
extern "C" const int const_GL_TEXTURE_DEPTH_SIZE;
extern "C" const int const_GL_TEXTURE_COMPARE_MODE;
extern "C" const int const_GL_TEXTURE_COMPARE_FUNC;
extern "C" const int const_GL_BUFFER_SIZE;
extern "C" const int const_GL_BUFFER_USAGE;
extern "C" const int const_GL_QUERY_COUNTER_BITS;
extern "C" const int const_GL_CURRENT_QUERY;
extern "C" const int const_GL_QUERY_RESULT;
extern "C" const int const_GL_QUERY_RESULT_AVAILABLE;
extern "C" const int const_GL_ARRAY_BUFFER;
extern "C" const int const_GL_ELEMENT_ARRAY_BUFFER;
extern "C" const int const_GL_ARRAY_BUFFER_BINDING;
extern "C" const int const_GL_ELEMENT_ARRAY_BUFFER_BINDING;
extern "C" const int const_GL_VERTEX_ATTRIB_ARRAY_BUFFER_BINDING;
extern "C" const int const_GL_READ_ONLY;
extern "C" const int const_GL_WRITE_ONLY;
extern "C" const int const_GL_READ_WRITE;
extern "C" const int const_GL_BUFFER_ACCESS;
extern "C" const int const_GL_BUFFER_MAPPED;
extern "C" const int const_GL_BUFFER_MAP_POINTER;
extern "C" const int const_GL_STREAM_DRAW;
extern "C" const int const_GL_STREAM_READ;
extern "C" const int const_GL_STREAM_COPY;
extern "C" const int const_GL_STATIC_DRAW;
extern "C" const int const_GL_STATIC_READ;
extern "C" const int const_GL_STATIC_COPY;
extern "C" const int const_GL_DYNAMIC_DRAW;
extern "C" const int const_GL_DYNAMIC_READ;
extern "C" const int const_GL_DYNAMIC_COPY;
extern "C" const int const_GL_SAMPLES_PASSED;
extern "C" const int const_GL_SRC1_ALPHA;
extern "C" const int const_GL_BLEND_EQUATION_RGB;
extern "C" const int const_GL_VERTEX_ATTRIB_ARRAY_ENABLED;
extern "C" const int const_GL_VERTEX_ATTRIB_ARRAY_SIZE;
extern "C" const int const_GL_VERTEX_ATTRIB_ARRAY_STRIDE;
extern "C" const int const_GL_VERTEX_ATTRIB_ARRAY_TYPE;
extern "C" const int const_GL_CURRENT_VERTEX_ATTRIB;
extern "C" const int const_GL_VERTEX_PROGRAM_POINT_SIZE;
extern "C" const int const_GL_VERTEX_ATTRIB_ARRAY_POINTER;
extern "C" const int const_GL_STENCIL_BACK_FUNC;
extern "C" const int const_GL_STENCIL_BACK_FAIL;
extern "C" const int const_GL_STENCIL_BACK_PASS_DEPTH_FAIL;
extern "C" const int const_GL_STENCIL_BACK_PASS_DEPTH_PASS;
extern "C" const int const_GL_MAX_DRAW_BUFFERS;
extern "C" const int const_GL_DRAW_BUFFER0;
extern "C" const int const_GL_DRAW_BUFFER1;
extern "C" const int const_GL_DRAW_BUFFER2;
extern "C" const int const_GL_DRAW_BUFFER3;
extern "C" const int const_GL_DRAW_BUFFER4;
extern "C" const int const_GL_DRAW_BUFFER5;
extern "C" const int const_GL_DRAW_BUFFER6;
extern "C" const int const_GL_DRAW_BUFFER7;
extern "C" const int const_GL_DRAW_BUFFER8;
extern "C" const int const_GL_DRAW_BUFFER9;
extern "C" const int const_GL_DRAW_BUFFER10;
extern "C" const int const_GL_DRAW_BUFFER11;
extern "C" const int const_GL_DRAW_BUFFER12;
extern "C" const int const_GL_DRAW_BUFFER13;
extern "C" const int const_GL_DRAW_BUFFER14;
extern "C" const int const_GL_DRAW_BUFFER15;
extern "C" const int const_GL_BLEND_EQUATION_ALPHA;
extern "C" const int const_GL_MAX_VERTEX_ATTRIBS;
extern "C" const int const_GL_VERTEX_ATTRIB_ARRAY_NORMALIZED;
extern "C" const int const_GL_MAX_TEXTURE_IMAGE_UNITS;
extern "C" const int const_GL_FRAGMENT_SHADER;
extern "C" const int const_GL_VERTEX_SHADER;
extern "C" const int const_GL_MAX_FRAGMENT_UNIFORM_COMPONENTS;
extern "C" const int const_GL_MAX_VERTEX_UNIFORM_COMPONENTS;
extern "C" const int const_GL_MAX_VARYING_FLOATS;
extern "C" const int const_GL_MAX_VERTEX_TEXTURE_IMAGE_UNITS;
extern "C" const int const_GL_MAX_COMBINED_TEXTURE_IMAGE_UNITS;
extern "C" const int const_GL_SHADER_TYPE;
extern "C" const int const_GL_FLOAT_VEC2;
extern "C" const int const_GL_FLOAT_VEC3;
extern "C" const int const_GL_FLOAT_VEC4;
extern "C" const int const_GL_INT_VEC2;
extern "C" const int const_GL_INT_VEC3;
extern "C" const int const_GL_INT_VEC4;
extern "C" const int const_GL_BOOL;
extern "C" const int const_GL_BOOL_VEC2;
extern "C" const int const_GL_BOOL_VEC3;
extern "C" const int const_GL_BOOL_VEC4;
extern "C" const int const_GL_FLOAT_MAT2;
extern "C" const int const_GL_FLOAT_MAT3;
extern "C" const int const_GL_FLOAT_MAT4;
extern "C" const int const_GL_SAMPLER_1D;
extern "C" const int const_GL_SAMPLER_2D;
extern "C" const int const_GL_SAMPLER_3D;
extern "C" const int const_GL_SAMPLER_CUBE;
extern "C" const int const_GL_SAMPLER_1D_SHADOW;
extern "C" const int const_GL_SAMPLER_2D_SHADOW;
extern "C" const int const_GL_DELETE_STATUS;
extern "C" const int const_GL_COMPILE_STATUS;
extern "C" const int const_GL_LINK_STATUS;
extern "C" const int const_GL_VALIDATE_STATUS;
extern "C" const int const_GL_INFO_LOG_LENGTH;
extern "C" const int const_GL_ATTACHED_SHADERS;
extern "C" const int const_GL_ACTIVE_UNIFORMS;
extern "C" const int const_GL_ACTIVE_UNIFORM_MAX_LENGTH;
extern "C" const int const_GL_SHADER_SOURCE_LENGTH;
extern "C" const int const_GL_ACTIVE_ATTRIBUTES;
extern "C" const int const_GL_ACTIVE_ATTRIBUTE_MAX_LENGTH;
extern "C" const int const_GL_FRAGMENT_SHADER_DERIVATIVE_HINT;
extern "C" const int const_GL_SHADING_LANGUAGE_VERSION;
extern "C" const int const_GL_CURRENT_PROGRAM;
extern "C" const int const_GL_POINT_SPRITE_COORD_ORIGIN;
extern "C" const int const_GL_LOWER_LEFT;
extern "C" const int const_GL_UPPER_LEFT;
extern "C" const int const_GL_STENCIL_BACK_REF;
extern "C" const int const_GL_STENCIL_BACK_VALUE_MASK;
extern "C" const int const_GL_STENCIL_BACK_WRITEMASK;
extern "C" const int const_GL_PIXEL_PACK_BUFFER;
extern "C" const int const_GL_PIXEL_UNPACK_BUFFER;
extern "C" const int const_GL_PIXEL_PACK_BUFFER_BINDING;
extern "C" const int const_GL_PIXEL_UNPACK_BUFFER_BINDING;
extern "C" const int const_GL_FLOAT_MAT2x3;
extern "C" const int const_GL_FLOAT_MAT2x4;
extern "C" const int const_GL_FLOAT_MAT3x2;
extern "C" const int const_GL_FLOAT_MAT3x4;
extern "C" const int const_GL_FLOAT_MAT4x2;
extern "C" const int const_GL_FLOAT_MAT4x3;
extern "C" const int const_GL_SRGB;
extern "C" const int const_GL_SRGB8;
extern "C" const int const_GL_SRGB_ALPHA;
extern "C" const int const_GL_SRGB8_ALPHA8;
extern "C" const int const_GL_COMPRESSED_SRGB;
extern "C" const int const_GL_COMPRESSED_SRGB_ALPHA;
extern "C" const int const_GL_COMPARE_REF_TO_TEXTURE;
extern "C" const int const_GL_CLIP_DISTANCE0;
extern "C" const int const_GL_CLIP_DISTANCE1;
extern "C" const int const_GL_CLIP_DISTANCE2;
extern "C" const int const_GL_CLIP_DISTANCE3;
extern "C" const int const_GL_CLIP_DISTANCE4;
extern "C" const int const_GL_CLIP_DISTANCE5;
extern "C" const int const_GL_CLIP_DISTANCE6;
extern "C" const int const_GL_CLIP_DISTANCE7;
extern "C" const int const_GL_MAX_CLIP_DISTANCES;
extern "C" const int const_GL_MAJOR_VERSION;
extern "C" const int const_GL_MINOR_VERSION;
extern "C" const int const_GL_NUM_EXTENSIONS;
extern "C" const int const_GL_CONTEXT_FLAGS;
extern "C" const int const_GL_COMPRESSED_RED;
extern "C" const int const_GL_COMPRESSED_RG;
extern "C" const int const_GL_CONTEXT_FLAG_FORWARD_COMPATIBLE_BIT;
extern "C" const int const_GL_RGBA32F;
extern "C" const int const_GL_RGB32F;
extern "C" const int const_GL_RGBA16F;
extern "C" const int const_GL_RGB16F;
extern "C" const int const_GL_VERTEX_ATTRIB_ARRAY_INTEGER;
extern "C" const int const_GL_MAX_ARRAY_TEXTURE_LAYERS;
extern "C" const int const_GL_MIN_PROGRAM_TEXEL_OFFSET;
extern "C" const int const_GL_MAX_PROGRAM_TEXEL_OFFSET;
extern "C" const int const_GL_CLAMP_READ_COLOR;
extern "C" const int const_GL_FIXED_ONLY;
extern "C" const int const_GL_MAX_VARYING_COMPONENTS;
extern "C" const int const_GL_TEXTURE_1D_ARRAY;
extern "C" const int const_GL_PROXY_TEXTURE_1D_ARRAY;
extern "C" const int const_GL_TEXTURE_2D_ARRAY;
extern "C" const int const_GL_PROXY_TEXTURE_2D_ARRAY;
extern "C" const int const_GL_TEXTURE_BINDING_1D_ARRAY;
extern "C" const int const_GL_TEXTURE_BINDING_2D_ARRAY;
extern "C" const int const_GL_R11F_G11F_B10F;
extern "C" const int const_GL_UNSIGNED_INT_10F_11F_11F_REV;
extern "C" const int const_GL_RGB9_E5;
extern "C" const int const_GL_UNSIGNED_INT_5_9_9_9_REV;
extern "C" const int const_GL_TEXTURE_SHARED_SIZE;
extern "C" const int const_GL_TRANSFORM_FEEDBACK_VARYING_MAX_LENGTH;
extern "C" const int const_GL_TRANSFORM_FEEDBACK_BUFFER_MODE;
extern "C" const int const_GL_MAX_TRANSFORM_FEEDBACK_SEPARATE_COMPONENTS;
extern "C" const int const_GL_TRANSFORM_FEEDBACK_VARYINGS;
extern "C" const int const_GL_TRANSFORM_FEEDBACK_BUFFER_START;
extern "C" const int const_GL_TRANSFORM_FEEDBACK_BUFFER_SIZE;
extern "C" const int const_GL_PRIMITIVES_GENERATED;
extern "C" const int const_GL_TRANSFORM_FEEDBACK_PRIMITIVES_WRITTEN;
extern "C" const int const_GL_RASTERIZER_DISCARD;
extern "C" const int const_GL_MAX_TRANSFORM_FEEDBACK_INTERLEAVED_COMPONENTS;
extern "C" const int const_GL_MAX_TRANSFORM_FEEDBACK_SEPARATE_ATTRIBS;
extern "C" const int const_GL_INTERLEAVED_ATTRIBS;
extern "C" const int const_GL_SEPARATE_ATTRIBS;
extern "C" const int const_GL_TRANSFORM_FEEDBACK_BUFFER;
extern "C" const int const_GL_TRANSFORM_FEEDBACK_BUFFER_BINDING;
extern "C" const int const_GL_RGBA32UI;
extern "C" const int const_GL_RGB32UI;
extern "C" const int const_GL_RGBA16UI;
extern "C" const int const_GL_RGB16UI;
extern "C" const int const_GL_RGBA8UI;
extern "C" const int const_GL_RGB8UI;
extern "C" const int const_GL_RGBA32I;
extern "C" const int const_GL_RGB32I;
extern "C" const int const_GL_RGBA16I;
extern "C" const int const_GL_RGB16I;
extern "C" const int const_GL_RGBA8I;
extern "C" const int const_GL_RGB8I;
extern "C" const int const_GL_RED_INTEGER;
extern "C" const int const_GL_GREEN_INTEGER;
extern "C" const int const_GL_BLUE_INTEGER;
extern "C" const int const_GL_RGB_INTEGER;
extern "C" const int const_GL_RGBA_INTEGER;
extern "C" const int const_GL_BGR_INTEGER;
extern "C" const int const_GL_BGRA_INTEGER;
extern "C" const int const_GL_SAMPLER_1D_ARRAY;
extern "C" const int const_GL_SAMPLER_2D_ARRAY;
extern "C" const int const_GL_SAMPLER_1D_ARRAY_SHADOW;
extern "C" const int const_GL_SAMPLER_2D_ARRAY_SHADOW;
extern "C" const int const_GL_SAMPLER_CUBE_SHADOW;
extern "C" const int const_GL_UNSIGNED_INT_VEC2;
extern "C" const int const_GL_UNSIGNED_INT_VEC3;
extern "C" const int const_GL_UNSIGNED_INT_VEC4;
extern "C" const int const_GL_INT_SAMPLER_1D;
extern "C" const int const_GL_INT_SAMPLER_2D;
extern "C" const int const_GL_INT_SAMPLER_3D;
extern "C" const int const_GL_INT_SAMPLER_CUBE;
extern "C" const int const_GL_INT_SAMPLER_1D_ARRAY;
extern "C" const int const_GL_INT_SAMPLER_2D_ARRAY;
extern "C" const int const_GL_UNSIGNED_INT_SAMPLER_1D;
extern "C" const int const_GL_UNSIGNED_INT_SAMPLER_2D;
extern "C" const int const_GL_UNSIGNED_INT_SAMPLER_3D;
extern "C" const int const_GL_UNSIGNED_INT_SAMPLER_CUBE;
extern "C" const int const_GL_UNSIGNED_INT_SAMPLER_1D_ARRAY;
extern "C" const int const_GL_UNSIGNED_INT_SAMPLER_2D_ARRAY;
extern "C" const int const_GL_QUERY_WAIT;
extern "C" const int const_GL_QUERY_NO_WAIT;
extern "C" const int const_GL_QUERY_BY_REGION_WAIT;
extern "C" const int const_GL_QUERY_BY_REGION_NO_WAIT;
extern "C" const int const_GL_BUFFER_ACCESS_FLAGS;
extern "C" const int const_GL_BUFFER_MAP_LENGTH;
extern "C" const int const_GL_BUFFER_MAP_OFFSET;
extern "C" const int const_GL_SAMPLER_2D_RECT;
extern "C" const int const_GL_SAMPLER_2D_RECT_SHADOW;
extern "C" const int const_GL_SAMPLER_BUFFER;
extern "C" const int const_GL_INT_SAMPLER_2D_RECT;
extern "C" const int const_GL_INT_SAMPLER_BUFFER;
extern "C" const int const_GL_UNSIGNED_INT_SAMPLER_2D_RECT;
extern "C" const int const_GL_UNSIGNED_INT_SAMPLER_BUFFER;
extern "C" const int const_GL_TEXTURE_BUFFER;
extern "C" const int const_GL_MAX_TEXTURE_BUFFER_SIZE;
extern "C" const int const_GL_TEXTURE_BINDING_BUFFER;
extern "C" const int const_GL_TEXTURE_BUFFER_DATA_STORE_BINDING;
extern "C" const int const_GL_TEXTURE_RECTANGLE;
extern "C" const int const_GL_TEXTURE_BINDING_RECTANGLE;
extern "C" const int const_GL_PROXY_TEXTURE_RECTANGLE;
extern "C" const int const_GL_MAX_RECTANGLE_TEXTURE_SIZE;
extern "C" const int const_GL_RED_SNORM;
extern "C" const int const_GL_RG_SNORM;
extern "C" const int const_GL_RGB_SNORM;
extern "C" const int const_GL_RGBA_SNORM;
extern "C" const int const_GL_R8_SNORM;
extern "C" const int const_GL_RG8_SNORM;
extern "C" const int const_GL_RGB8_SNORM;
extern "C" const int const_GL_RGBA8_SNORM;
extern "C" const int const_GL_R16_SNORM;
extern "C" const int const_GL_RG16_SNORM;
extern "C" const int const_GL_RGB16_SNORM;
extern "C" const int const_GL_RGBA16_SNORM;
extern "C" const int const_GL_SIGNED_NORMALIZED;
extern "C" const int const_GL_PRIMITIVE_RESTART;
extern "C" const int const_GL_PRIMITIVE_RESTART_INDEX;
extern "C" const int const_GL_CONTEXT_CORE_PROFILE_BIT;
extern "C" const int const_GL_CONTEXT_COMPATIBILITY_PROFILE_BIT;
extern "C" const int const_GL_LINES_ADJACENCY;
extern "C" const int const_GL_LINE_STRIP_ADJACENCY;
extern "C" const int const_GL_TRIANGLES_ADJACENCY;
extern "C" const int const_GL_TRIANGLE_STRIP_ADJACENCY;
extern "C" const int const_GL_PROGRAM_POINT_SIZE;
extern "C" const int const_GL_MAX_GEOMETRY_TEXTURE_IMAGE_UNITS;
extern "C" const int const_GL_FRAMEBUFFER_ATTACHMENT_LAYERED;
extern "C" const int const_GL_FRAMEBUFFER_INCOMPLETE_LAYER_TARGETS;
extern "C" const int const_GL_GEOMETRY_SHADER;
extern "C" const int const_GL_GEOMETRY_VERTICES_OUT;
extern "C" const int const_GL_GEOMETRY_INPUT_TYPE;
extern "C" const int const_GL_GEOMETRY_OUTPUT_TYPE;
extern "C" const int const_GL_MAX_GEOMETRY_UNIFORM_COMPONENTS;
extern "C" const int const_GL_MAX_GEOMETRY_OUTPUT_VERTICES;
extern "C" const int const_GL_MAX_GEOMETRY_TOTAL_OUTPUT_COMPONENTS;
extern "C" const int const_GL_MAX_VERTEX_OUTPUT_COMPONENTS;
extern "C" const int const_GL_MAX_GEOMETRY_INPUT_COMPONENTS;
extern "C" const int const_GL_MAX_GEOMETRY_OUTPUT_COMPONENTS;
extern "C" const int const_GL_MAX_FRAGMENT_INPUT_COMPONENTS;
extern "C" const int const_GL_CONTEXT_PROFILE_MASK;
extern "C" const int const_GL_VERTEX_ATTRIB_ARRAY_DIVISOR;
extern "C" const int const_GL_SAMPLE_SHADING;
extern "C" const int const_GL_MIN_SAMPLE_SHADING_VALUE;
extern "C" const int const_GL_MIN_PROGRAM_TEXTURE_GATHER_OFFSET;
extern "C" const int const_GL_MAX_PROGRAM_TEXTURE_GATHER_OFFSET;
extern "C" const int const_GL_TEXTURE_CUBE_MAP_ARRAY;
extern "C" const int const_GL_TEXTURE_BINDING_CUBE_MAP_ARRAY;
extern "C" const int const_GL_PROXY_TEXTURE_CUBE_MAP_ARRAY;
extern "C" const int const_GL_SAMPLER_CUBE_MAP_ARRAY;
extern "C" const int const_GL_SAMPLER_CUBE_MAP_ARRAY_SHADOW;
extern "C" const int const_GL_INT_SAMPLER_CUBE_MAP_ARRAY;
extern "C" const int const_GL_UNSIGNED_INT_SAMPLER_CUBE_MAP_ARRAY;
extern "C" const int const_GL_NUM_SHADING_LANGUAGE_VERSIONS;
extern "C" const int const_GL_VERTEX_ATTRIB_ARRAY_LONG;
extern "C" const int const_GL_DEPTH_COMPONENT32F;
extern "C" const int const_GL_DEPTH32F_STENCIL8;
extern "C" const int const_GL_FLOAT_32_UNSIGNED_INT_24_8_REV;
extern "C" const int const_GL_INVALID_FRAMEBUFFER_OPERATION;
extern "C" const int const_GL_FRAMEBUFFER_ATTACHMENT_COLOR_ENCODING;
extern "C" const int const_GL_FRAMEBUFFER_ATTACHMENT_COMPONENT_TYPE;
extern "C" const int const_GL_FRAMEBUFFER_ATTACHMENT_RED_SIZE;
extern "C" const int const_GL_FRAMEBUFFER_ATTACHMENT_GREEN_SIZE;
extern "C" const int const_GL_FRAMEBUFFER_ATTACHMENT_BLUE_SIZE;
extern "C" const int const_GL_FRAMEBUFFER_ATTACHMENT_ALPHA_SIZE;
extern "C" const int const_GL_FRAMEBUFFER_ATTACHMENT_DEPTH_SIZE;
extern "C" const int const_GL_FRAMEBUFFER_ATTACHMENT_STENCIL_SIZE;
extern "C" const int const_GL_FRAMEBUFFER_DEFAULT;
extern "C" const int const_GL_FRAMEBUFFER_UNDEFINED;
extern "C" const int const_GL_DEPTH_STENCIL_ATTACHMENT;
extern "C" const int const_GL_MAX_RENDERBUFFER_SIZE;
extern "C" const int const_GL_DEPTH_STENCIL;
extern "C" const int const_GL_UNSIGNED_INT_24_8;
extern "C" const int const_GL_DEPTH24_STENCIL8;
extern "C" const int const_GL_TEXTURE_STENCIL_SIZE;
extern "C" const int const_GL_TEXTURE_RED_TYPE;
extern "C" const int const_GL_TEXTURE_GREEN_TYPE;
extern "C" const int const_GL_TEXTURE_BLUE_TYPE;
extern "C" const int const_GL_TEXTURE_ALPHA_TYPE;
extern "C" const int const_GL_TEXTURE_DEPTH_TYPE;
extern "C" const int const_GL_UNSIGNED_NORMALIZED;
extern "C" const int const_GL_FRAMEBUFFER_BINDING;
extern "C" const int const_GL_DRAW_FRAMEBUFFER_BINDING;
extern "C" const int const_GL_RENDERBUFFER_BINDING;
extern "C" const int const_GL_READ_FRAMEBUFFER;
extern "C" const int const_GL_DRAW_FRAMEBUFFER;
extern "C" const int const_GL_READ_FRAMEBUFFER_BINDING;
extern "C" const int const_GL_RENDERBUFFER_SAMPLES;
extern "C" const int const_GL_FRAMEBUFFER_ATTACHMENT_OBJECT_TYPE;
extern "C" const int const_GL_FRAMEBUFFER_ATTACHMENT_OBJECT_NAME;
extern "C" const int const_GL_FRAMEBUFFER_ATTACHMENT_TEXTURE_LEVEL;
extern "C" const int const_GL_FRAMEBUFFER_ATTACHMENT_TEXTURE_CUBE_MAP_FACE;
extern "C" const int const_GL_FRAMEBUFFER_ATTACHMENT_TEXTURE_LAYER;
extern "C" const int const_GL_FRAMEBUFFER_COMPLETE;
extern "C" const int const_GL_FRAMEBUFFER_INCOMPLETE_ATTACHMENT;
extern "C" const int const_GL_FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT;
extern "C" const int const_GL_FRAMEBUFFER_INCOMPLETE_DRAW_BUFFER;
extern "C" const int const_GL_FRAMEBUFFER_INCOMPLETE_READ_BUFFER;
extern "C" const int const_GL_FRAMEBUFFER_UNSUPPORTED;
extern "C" const int const_GL_MAX_COLOR_ATTACHMENTS;
extern "C" const int const_GL_COLOR_ATTACHMENT0;
extern "C" const int const_GL_COLOR_ATTACHMENT1;
extern "C" const int const_GL_COLOR_ATTACHMENT2;
extern "C" const int const_GL_COLOR_ATTACHMENT3;
extern "C" const int const_GL_COLOR_ATTACHMENT4;
extern "C" const int const_GL_COLOR_ATTACHMENT5;
extern "C" const int const_GL_COLOR_ATTACHMENT6;
extern "C" const int const_GL_COLOR_ATTACHMENT7;
extern "C" const int const_GL_COLOR_ATTACHMENT8;
extern "C" const int const_GL_COLOR_ATTACHMENT9;
extern "C" const int const_GL_COLOR_ATTACHMENT10;
extern "C" const int const_GL_COLOR_ATTACHMENT11;
extern "C" const int const_GL_COLOR_ATTACHMENT12;
extern "C" const int const_GL_COLOR_ATTACHMENT13;
extern "C" const int const_GL_COLOR_ATTACHMENT14;
extern "C" const int const_GL_COLOR_ATTACHMENT15;
extern "C" const int const_GL_DEPTH_ATTACHMENT;
extern "C" const int const_GL_STENCIL_ATTACHMENT;
extern "C" const int const_GL_FRAMEBUFFER;
extern "C" const int const_GL_RENDERBUFFER;
extern "C" const int const_GL_RENDERBUFFER_WIDTH;
extern "C" const int const_GL_RENDERBUFFER_HEIGHT;
extern "C" const int const_GL_RENDERBUFFER_INTERNAL_FORMAT;
extern "C" const int const_GL_STENCIL_INDEX1;
extern "C" const int const_GL_STENCIL_INDEX4;
extern "C" const int const_GL_STENCIL_INDEX8;
extern "C" const int const_GL_STENCIL_INDEX16;
extern "C" const int const_GL_RENDERBUFFER_RED_SIZE;
extern "C" const int const_GL_RENDERBUFFER_GREEN_SIZE;
extern "C" const int const_GL_RENDERBUFFER_BLUE_SIZE;
extern "C" const int const_GL_RENDERBUFFER_ALPHA_SIZE;
extern "C" const int const_GL_RENDERBUFFER_DEPTH_SIZE;
extern "C" const int const_GL_RENDERBUFFER_STENCIL_SIZE;
extern "C" const int const_GL_FRAMEBUFFER_INCOMPLETE_MULTISAMPLE;
extern "C" const int const_GL_MAX_SAMPLES;
extern "C" const int const_GL_FRAMEBUFFER_SRGB;
extern "C" const int const_GL_HALF_FLOAT;
extern "C" const int const_GL_MAP_READ_BIT;
extern "C" const int const_GL_MAP_WRITE_BIT;
extern "C" const int const_GL_MAP_INVALIDATE_RANGE_BIT;
extern "C" const int const_GL_MAP_INVALIDATE_BUFFER_BIT;
extern "C" const int const_GL_MAP_FLUSH_EXPLICIT_BIT;
extern "C" const int const_GL_MAP_UNSYNCHRONIZED_BIT;
extern "C" const int const_GL_COMPRESSED_RED_RGTC1;
extern "C" const int const_GL_COMPRESSED_SIGNED_RED_RGTC1;
extern "C" const int const_GL_COMPRESSED_RG_RGTC2;
extern "C" const int const_GL_COMPRESSED_SIGNED_RG_RGTC2;
extern "C" const int const_GL_RG;
extern "C" const int const_GL_RG_INTEGER;
extern "C" const int const_GL_R8;
extern "C" const int const_GL_R16;
extern "C" const int const_GL_RG8;
extern "C" const int const_GL_RG16;
extern "C" const int const_GL_R16F;
extern "C" const int const_GL_R32F;
extern "C" const int const_GL_RG16F;
extern "C" const int const_GL_RG32F;
extern "C" const int const_GL_R8I;
extern "C" const int const_GL_R8UI;
extern "C" const int const_GL_R16I;
extern "C" const int const_GL_R16UI;
extern "C" const int const_GL_R32I;
extern "C" const int const_GL_R32UI;
extern "C" const int const_GL_RG8I;
extern "C" const int const_GL_RG8UI;
extern "C" const int const_GL_RG16I;
extern "C" const int const_GL_RG16UI;
extern "C" const int const_GL_RG32I;
extern "C" const int const_GL_RG32UI;
extern "C" const int const_GL_VERTEX_ARRAY_BINDING;
extern "C" const int const_GL_UNIFORM_BUFFER;
extern "C" const int const_GL_UNIFORM_BUFFER_BINDING;
extern "C" const int const_GL_UNIFORM_BUFFER_START;
extern "C" const int const_GL_UNIFORM_BUFFER_SIZE;
extern "C" const int const_GL_MAX_VERTEX_UNIFORM_BLOCKS;
extern "C" const int const_GL_MAX_GEOMETRY_UNIFORM_BLOCKS;
extern "C" const int const_GL_MAX_FRAGMENT_UNIFORM_BLOCKS;
extern "C" const int const_GL_MAX_COMBINED_UNIFORM_BLOCKS;
extern "C" const int const_GL_MAX_UNIFORM_BUFFER_BINDINGS;
extern "C" const int const_GL_MAX_UNIFORM_BLOCK_SIZE;
extern "C" const int const_GL_MAX_COMBINED_VERTEX_UNIFORM_COMPONENTS;
extern "C" const int const_GL_MAX_COMBINED_GEOMETRY_UNIFORM_COMPONENTS;
extern "C" const int const_GL_MAX_COMBINED_FRAGMENT_UNIFORM_COMPONENTS;
extern "C" const int const_GL_UNIFORM_BUFFER_OFFSET_ALIGNMENT;
extern "C" const int const_GL_ACTIVE_UNIFORM_BLOCK_MAX_NAME_LENGTH;
extern "C" const int const_GL_ACTIVE_UNIFORM_BLOCKS;
extern "C" const int const_GL_UNIFORM_TYPE;
extern "C" const int const_GL_UNIFORM_SIZE;
extern "C" const int const_GL_UNIFORM_NAME_LENGTH;
extern "C" const int const_GL_UNIFORM_BLOCK_INDEX;
extern "C" const int const_GL_UNIFORM_OFFSET;
extern "C" const int const_GL_UNIFORM_ARRAY_STRIDE;
extern "C" const int const_GL_UNIFORM_MATRIX_STRIDE;
extern "C" const int const_GL_UNIFORM_IS_ROW_MAJOR;
extern "C" const int const_GL_UNIFORM_BLOCK_BINDING;
extern "C" const int const_GL_UNIFORM_BLOCK_DATA_SIZE;
extern "C" const int const_GL_UNIFORM_BLOCK_NAME_LENGTH;
extern "C" const int const_GL_UNIFORM_BLOCK_ACTIVE_UNIFORMS;
extern "C" const int const_GL_UNIFORM_BLOCK_ACTIVE_UNIFORM_INDICES;
extern "C" const int const_GL_UNIFORM_BLOCK_REFERENCED_BY_VERTEX_SHADER;
extern "C" const int const_GL_UNIFORM_BLOCK_REFERENCED_BY_GEOMETRY_SHADER;
extern "C" const int const_GL_UNIFORM_BLOCK_REFERENCED_BY_FRAGMENT_SHADER;
extern "C" const int const_GL_INVALID_INDEX;
extern "C" const int const_GL_COPY_READ_BUFFER_BINDING;
extern "C" const int const_GL_COPY_READ_BUFFER;
extern "C" const int const_GL_COPY_WRITE_BUFFER_BINDING;
extern "C" const int const_GL_COPY_WRITE_BUFFER;
extern "C" const int const_GL_DEPTH_CLAMP;
extern "C" const int const_GL_QUADS_FOLLOW_PROVOKING_VERTEX_CONVENTION;
extern "C" const int const_GL_FIRST_VERTEX_CONVENTION;
extern "C" const int const_GL_LAST_VERTEX_CONVENTION;
extern "C" const int const_GL_PROVOKING_VERTEX;
extern "C" const int const_GL_TEXTURE_CUBE_MAP_SEAMLESS;
extern "C" const int const_GL_MAX_SERVER_WAIT_TIMEOUT;
extern "C" const int const_GL_OBJECT_TYPE;
extern "C" const int const_GL_SYNC_CONDITION;
extern "C" const int const_GL_SYNC_STATUS;
extern "C" const int const_GL_SYNC_FLAGS;
extern "C" const int const_GL_SYNC_FENCE;
extern "C" const int const_GL_SYNC_GPU_COMMANDS_COMPLETE;
extern "C" const int const_GL_UNSIGNALED;
extern "C" const int const_GL_SIGNALED;
extern "C" const int const_GL_ALREADY_SIGNALED;
extern "C" const int const_GL_TIMEOUT_EXPIRED;
extern "C" const int const_GL_CONDITION_SATISFIED;
extern "C" const int const_GL_WAIT_FAILED;
extern "C" const int const_GL_SYNC_FLUSH_COMMANDS_BIT;
extern "C" const int const_GL_TIMEOUT_IGNORED;
extern "C" const int const_GL_SAMPLE_POSITION;
extern "C" const int const_GL_SAMPLE_MASK;
extern "C" const int const_GL_SAMPLE_MASK_VALUE;
extern "C" const int const_GL_MAX_SAMPLE_MASK_WORDS;
extern "C" const int const_GL_TEXTURE_2D_MULTISAMPLE;
extern "C" const int const_GL_PROXY_TEXTURE_2D_MULTISAMPLE;
extern "C" const int const_GL_TEXTURE_2D_MULTISAMPLE_ARRAY;
extern "C" const int const_GL_PROXY_TEXTURE_2D_MULTISAMPLE_ARRAY;
extern "C" const int const_GL_TEXTURE_BINDING_2D_MULTISAMPLE;
extern "C" const int const_GL_TEXTURE_BINDING_2D_MULTISAMPLE_ARRAY;
extern "C" const int const_GL_TEXTURE_SAMPLES;
extern "C" const int const_GL_TEXTURE_FIXED_SAMPLE_LOCATIONS;
extern "C" const int const_GL_SAMPLER_2D_MULTISAMPLE;
extern "C" const int const_GL_INT_SAMPLER_2D_MULTISAMPLE;
extern "C" const int const_GL_UNSIGNED_INT_SAMPLER_2D_MULTISAMPLE;
extern "C" const int const_GL_SAMPLER_2D_MULTISAMPLE_ARRAY;
extern "C" const int const_GL_INT_SAMPLER_2D_MULTISAMPLE_ARRAY;
extern "C" const int const_GL_UNSIGNED_INT_SAMPLER_2D_MULTISAMPLE_ARRAY;
extern "C" const int const_GL_MAX_COLOR_TEXTURE_SAMPLES;
extern "C" const int const_GL_MAX_DEPTH_TEXTURE_SAMPLES;
extern "C" const int const_GL_MAX_INTEGER_SAMPLES;
extern "C" const int const_GL_SAMPLE_SHADING_ARB;
extern "C" const int const_GL_MIN_SAMPLE_SHADING_VALUE_ARB;
extern "C" const int const_GL_TEXTURE_CUBE_MAP_ARRAY_ARB;
extern "C" const int const_GL_TEXTURE_BINDING_CUBE_MAP_ARRAY_ARB;
extern "C" const int const_GL_PROXY_TEXTURE_CUBE_MAP_ARRAY_ARB;
extern "C" const int const_GL_SAMPLER_CUBE_MAP_ARRAY_ARB;
extern "C" const int const_GL_SAMPLER_CUBE_MAP_ARRAY_SHADOW_ARB;
extern "C" const int const_GL_INT_SAMPLER_CUBE_MAP_ARRAY_ARB;
extern "C" const int const_GL_UNSIGNED_INT_SAMPLER_CUBE_MAP_ARRAY_ARB;
extern "C" const int const_GL_MIN_PROGRAM_TEXTURE_GATHER_OFFSET_ARB;
extern "C" const int const_GL_MAX_PROGRAM_TEXTURE_GATHER_OFFSET_ARB;
extern "C" const int const_GL_MAX_PROGRAM_TEXTURE_GATHER_COMPONENTS_ARB;
extern "C" const int const_GL_SHADER_INCLUDE_ARB;
extern "C" const int const_GL_NAMED_STRING_LENGTH_ARB;
extern "C" const int const_GL_NAMED_STRING_TYPE_ARB;
extern "C" const int const_GL_COMPRESSED_RGBA_BPTC_UNORM_ARB;
extern "C" const int const_GL_COMPRESSED_SRGB_ALPHA_BPTC_UNORM_ARB;
extern "C" const int const_GL_COMPRESSED_RGB_BPTC_SIGNED_FLOAT_ARB;
extern "C" const int const_GL_COMPRESSED_RGB_BPTC_UNSIGNED_FLOAT_ARB;
extern "C" const int const_GL_SRC1_COLOR;
extern "C" const int const_GL_ONE_MINUS_SRC1_COLOR;
extern "C" const int const_GL_ONE_MINUS_SRC1_ALPHA;
extern "C" const int const_GL_MAX_DUAL_SOURCE_DRAW_BUFFERS;
extern "C" const int const_GL_ANY_SAMPLES_PASSED;
extern "C" const int const_GL_SAMPLER_BINDING;
extern "C" const int const_GL_RGB10_A2UI;
extern "C" const int const_GL_TEXTURE_SWIZZLE_R;
extern "C" const int const_GL_TEXTURE_SWIZZLE_G;
extern "C" const int const_GL_TEXTURE_SWIZZLE_B;
extern "C" const int const_GL_TEXTURE_SWIZZLE_A;
extern "C" const int const_GL_TEXTURE_SWIZZLE_RGBA;
extern "C" const int const_GL_TIME_ELAPSED;
extern "C" const int const_GL_TIMESTAMP;
extern "C" const int const_GL_INT_2_10_10_10_REV;
extern "C" const int const_GL_DRAW_INDIRECT_BUFFER;
extern "C" const int const_GL_DRAW_INDIRECT_BUFFER_BINDING;
extern "C" const int const_GL_GEOMETRY_SHADER_INVOCATIONS;
extern "C" const int const_GL_MAX_GEOMETRY_SHADER_INVOCATIONS;
extern "C" const int const_GL_MIN_FRAGMENT_INTERPOLATION_OFFSET;
extern "C" const int const_GL_MAX_FRAGMENT_INTERPOLATION_OFFSET;
extern "C" const int const_GL_FRAGMENT_INTERPOLATION_OFFSET_BITS;
extern "C" const int const_GL_DOUBLE_VEC2;
extern "C" const int const_GL_DOUBLE_VEC3;
extern "C" const int const_GL_DOUBLE_VEC4;
extern "C" const int const_GL_DOUBLE_MAT2;
extern "C" const int const_GL_DOUBLE_MAT3;
extern "C" const int const_GL_DOUBLE_MAT4;
extern "C" const int const_GL_DOUBLE_MAT2x3;
extern "C" const int const_GL_DOUBLE_MAT2x4;
extern "C" const int const_GL_DOUBLE_MAT3x2;
extern "C" const int const_GL_DOUBLE_MAT3x4;
extern "C" const int const_GL_DOUBLE_MAT4x2;
extern "C" const int const_GL_DOUBLE_MAT4x3;
extern "C" const int const_GL_ACTIVE_SUBROUTINES;
extern "C" const int const_GL_ACTIVE_SUBROUTINE_UNIFORMS;
extern "C" const int const_GL_ACTIVE_SUBROUTINE_UNIFORM_LOCATIONS;
extern "C" const int const_GL_ACTIVE_SUBROUTINE_MAX_LENGTH;
extern "C" const int const_GL_ACTIVE_SUBROUTINE_UNIFORM_MAX_LENGTH;
extern "C" const int const_GL_MAX_SUBROUTINES;
extern "C" const int const_GL_MAX_SUBROUTINE_UNIFORM_LOCATIONS;
extern "C" const int const_GL_NUM_COMPATIBLE_SUBROUTINES;
extern "C" const int const_GL_COMPATIBLE_SUBROUTINES;
extern "C" const int const_GL_PATCHES;
extern "C" const int const_GL_PATCH_VERTICES;
extern "C" const int const_GL_PATCH_DEFAULT_INNER_LEVEL;
extern "C" const int const_GL_PATCH_DEFAULT_OUTER_LEVEL;
extern "C" const int const_GL_TESS_CONTROL_OUTPUT_VERTICES;
extern "C" const int const_GL_TESS_GEN_MODE;
extern "C" const int const_GL_TESS_GEN_SPACING;
extern "C" const int const_GL_TESS_GEN_VERTEX_ORDER;
extern "C" const int const_GL_TESS_GEN_POINT_MODE;
extern "C" const int const_GL_ISOLINES;
extern "C" const int const_GL_FRACTIONAL_ODD;
extern "C" const int const_GL_FRACTIONAL_EVEN;
extern "C" const int const_GL_MAX_PATCH_VERTICES;
extern "C" const int const_GL_MAX_TESS_GEN_LEVEL;
extern "C" const int const_GL_MAX_TESS_CONTROL_UNIFORM_COMPONENTS;
extern "C" const int const_GL_MAX_TESS_EVALUATION_UNIFORM_COMPONENTS;
extern "C" const int const_GL_MAX_TESS_CONTROL_TEXTURE_IMAGE_UNITS;
extern "C" const int const_GL_MAX_TESS_EVALUATION_TEXTURE_IMAGE_UNITS;
extern "C" const int const_GL_MAX_TESS_CONTROL_OUTPUT_COMPONENTS;
extern "C" const int const_GL_MAX_TESS_PATCH_COMPONENTS;
extern "C" const int const_GL_MAX_TESS_CONTROL_TOTAL_OUTPUT_COMPONENTS;
extern "C" const int const_GL_MAX_TESS_EVALUATION_OUTPUT_COMPONENTS;
extern "C" const int const_GL_MAX_TESS_CONTROL_UNIFORM_BLOCKS;
extern "C" const int const_GL_MAX_TESS_EVALUATION_UNIFORM_BLOCKS;
extern "C" const int const_GL_MAX_TESS_CONTROL_INPUT_COMPONENTS;
extern "C" const int const_GL_MAX_TESS_EVALUATION_INPUT_COMPONENTS;
extern "C" const int const_GL_MAX_COMBINED_TESS_CONTROL_UNIFORM_COMPONENTS;
extern "C" const int const_GL_MAX_COMBINED_TESS_EVALUATION_UNIFORM_COMPONENTS;
extern "C" const int const_GL_UNIFORM_BLOCK_REFERENCED_BY_TESS_CONTROL_SHADER;
extern "C" const int const_GL_UNIFORM_BLOCK_REFERENCED_BY_TESS_EVALUATION_SHADER;
extern "C" const int const_GL_TESS_EVALUATION_SHADER;
extern "C" const int const_GL_TESS_CONTROL_SHADER;
extern "C" const int const_GL_TRANSFORM_FEEDBACK;
extern "C" const int const_GL_TRANSFORM_FEEDBACK_PAUSED;
extern "C" const int const_GL_TRANSFORM_FEEDBACK_BUFFER_PAUSED;
extern "C" const int const_GL_TRANSFORM_FEEDBACK_ACTIVE;
extern "C" const int const_GL_TRANSFORM_FEEDBACK_BUFFER_ACTIVE;
extern "C" const int const_GL_TRANSFORM_FEEDBACK_BINDING;
extern "C" const int const_GL_MAX_TRANSFORM_FEEDBACK_BUFFERS;
extern "C" const int const_GL_MAX_VERTEX_STREAMS;
extern "C" const int const_GL_FIXED;
extern "C" const int const_GL_IMPLEMENTATION_COLOR_READ_TYPE;
extern "C" const int const_GL_IMPLEMENTATION_COLOR_READ_FORMAT;
extern "C" const int const_GL_LOW_FLOAT;
extern "C" const int const_GL_MEDIUM_FLOAT;
extern "C" const int const_GL_HIGH_FLOAT;
extern "C" const int const_GL_LOW_INT;
extern "C" const int const_GL_MEDIUM_INT;
extern "C" const int const_GL_HIGH_INT;
extern "C" const int const_GL_SHADER_COMPILER;
extern "C" const int const_GL_SHADER_BINARY_FORMATS;
extern "C" const int const_GL_NUM_SHADER_BINARY_FORMATS;
extern "C" const int const_GL_MAX_VERTEX_UNIFORM_VECTORS;
extern "C" const int const_GL_MAX_VARYING_VECTORS;
extern "C" const int const_GL_MAX_FRAGMENT_UNIFORM_VECTORS;
extern "C" const int const_GL_RGB565;
extern "C" const int const_GL_PROGRAM_BINARY_RETRIEVABLE_HINT;
extern "C" const int const_GL_PROGRAM_BINARY_LENGTH;
extern "C" const int const_GL_NUM_PROGRAM_BINARY_FORMATS;
extern "C" const int const_GL_PROGRAM_BINARY_FORMATS;
extern "C" const int const_GL_VERTEX_SHADER_BIT;
extern "C" const int const_GL_FRAGMENT_SHADER_BIT;
extern "C" const int const_GL_GEOMETRY_SHADER_BIT;
extern "C" const int const_GL_TESS_CONTROL_SHADER_BIT;
extern "C" const int const_GL_TESS_EVALUATION_SHADER_BIT;
extern "C" const int const_GL_ALL_SHADER_BITS;
extern "C" const int const_GL_PROGRAM_SEPARABLE;
extern "C" const int const_GL_ACTIVE_PROGRAM;
extern "C" const int const_GL_PROGRAM_PIPELINE_BINDING;
extern "C" const int const_GL_MAX_VIEWPORTS;
extern "C" const int const_GL_VIEWPORT_SUBPIXEL_BITS;
extern "C" const int const_GL_VIEWPORT_BOUNDS_RANGE;
extern "C" const int const_GL_LAYER_PROVOKING_VERTEX;
extern "C" const int const_GL_VIEWPORT_INDEX_PROVOKING_VERTEX;
extern "C" const int const_GL_UNDEFINED_VERTEX;
extern "C" const int const_GL_SYNC_CL_EVENT_ARB;
extern "C" const int const_GL_SYNC_CL_EVENT_COMPLETE_ARB;
extern "C" const int const_GL_DEBUG_OUTPUT_SYNCHRONOUS_ARB;
extern "C" const int const_GL_DEBUG_NEXT_LOGGED_MESSAGE_LENGTH_ARB;
extern "C" const int const_GL_DEBUG_CALLBACK_FUNCTION_ARB;
extern "C" const int const_GL_DEBUG_CALLBACK_USER_PARAM_ARB;
extern "C" const int const_GL_DEBUG_SOURCE_API_ARB;
extern "C" const int const_GL_DEBUG_SOURCE_WINDOW_SYSTEM_ARB;
extern "C" const int const_GL_DEBUG_SOURCE_SHADER_COMPILER_ARB;
extern "C" const int const_GL_DEBUG_SOURCE_THIRD_PARTY_ARB;
extern "C" const int const_GL_DEBUG_SOURCE_APPLICATION_ARB;
extern "C" const int const_GL_DEBUG_SOURCE_OTHER_ARB;
extern "C" const int const_GL_DEBUG_TYPE_ERROR_ARB;
extern "C" const int const_GL_DEBUG_TYPE_DEPRECATED_BEHAVIOR_ARB;
extern "C" const int const_GL_DEBUG_TYPE_UNDEFINED_BEHAVIOR_ARB;
extern "C" const int const_GL_DEBUG_TYPE_PORTABILITY_ARB;
extern "C" const int const_GL_DEBUG_TYPE_PERFORMANCE_ARB;
extern "C" const int const_GL_DEBUG_TYPE_OTHER_ARB;
extern "C" const int const_GL_MAX_DEBUG_MESSAGE_LENGTH_ARB;
extern "C" const int const_GL_MAX_DEBUG_LOGGED_MESSAGES_ARB;
extern "C" const int const_GL_DEBUG_LOGGED_MESSAGES_ARB;
extern "C" const int const_GL_DEBUG_SEVERITY_HIGH_ARB;
extern "C" const int const_GL_DEBUG_SEVERITY_MEDIUM_ARB;
extern "C" const int const_GL_DEBUG_SEVERITY_LOW_ARB;
extern "C" const int const_GL_CONTEXT_FLAG_ROBUST_ACCESS_BIT_ARB;
extern "C" const int const_GL_LOSE_CONTEXT_ON_RESET_ARB;
extern "C" const int const_GL_GUILTY_CONTEXT_RESET_ARB;
extern "C" const int const_GL_INNOCENT_CONTEXT_RESET_ARB;
extern "C" const int const_GL_UNKNOWN_CONTEXT_RESET_ARB;
extern "C" const int const_GL_RESET_NOTIFICATION_STRATEGY_ARB;
extern "C" const int const_GL_NO_RESET_NOTIFICATION_ARB;
extern "C" const int const_GL_UNPACK_COMPRESSED_BLOCK_WIDTH;
extern "C" const int const_GL_UNPACK_COMPRESSED_BLOCK_HEIGHT;
extern "C" const int const_GL_UNPACK_COMPRESSED_BLOCK_DEPTH;
extern "C" const int const_GL_UNPACK_COMPRESSED_BLOCK_SIZE;
extern "C" const int const_GL_PACK_COMPRESSED_BLOCK_WIDTH;
extern "C" const int const_GL_PACK_COMPRESSED_BLOCK_HEIGHT;
extern "C" const int const_GL_PACK_COMPRESSED_BLOCK_DEPTH;
extern "C" const int const_GL_PACK_COMPRESSED_BLOCK_SIZE;
extern "C" const int const_GL_NUM_SAMPLE_COUNTS;
extern "C" const int const_GL_MIN_MAP_BUFFER_ALIGNMENT;
extern "C" const int const_GL_ATOMIC_COUNTER_BUFFER;
extern "C" const int const_GL_ATOMIC_COUNTER_BUFFER_BINDING;
extern "C" const int const_GL_ATOMIC_COUNTER_BUFFER_START;
extern "C" const int const_GL_ATOMIC_COUNTER_BUFFER_SIZE;
extern "C" const int const_GL_ATOMIC_COUNTER_BUFFER_DATA_SIZE;
extern "C" const int const_GL_ATOMIC_COUNTER_BUFFER_ACTIVE_ATOMIC_COUNTERS;
extern "C" const int const_GL_ATOMIC_COUNTER_BUFFER_ACTIVE_ATOMIC_COUNTER_INDICES;
extern "C" const int const_GL_ATOMIC_COUNTER_BUFFER_REFERENCED_BY_VERTEX_SHADER;
extern "C" const int const_GL_ATOMIC_COUNTER_BUFFER_REFERENCED_BY_TESS_CONTROL_SHADER;
extern "C" const int const_GL_ATOMIC_COUNTER_BUFFER_REFERENCED_BY_TESS_EVALUATION_SHADER;
extern "C" const int const_GL_ATOMIC_COUNTER_BUFFER_REFERENCED_BY_GEOMETRY_SHADER;
extern "C" const int const_GL_ATOMIC_COUNTER_BUFFER_REFERENCED_BY_FRAGMENT_SHADER;
extern "C" const int const_GL_MAX_VERTEX_ATOMIC_COUNTER_BUFFERS;
extern "C" const int const_GL_MAX_TESS_CONTROL_ATOMIC_COUNTER_BUFFERS;
extern "C" const int const_GL_MAX_TESS_EVALUATION_ATOMIC_COUNTER_BUFFERS;
extern "C" const int const_GL_MAX_GEOMETRY_ATOMIC_COUNTER_BUFFERS;
extern "C" const int const_GL_MAX_FRAGMENT_ATOMIC_COUNTER_BUFFERS;
extern "C" const int const_GL_MAX_COMBINED_ATOMIC_COUNTER_BUFFERS;
extern "C" const int const_GL_MAX_VERTEX_ATOMIC_COUNTERS;
extern "C" const int const_GL_MAX_TESS_CONTROL_ATOMIC_COUNTERS;
extern "C" const int const_GL_MAX_TESS_EVALUATION_ATOMIC_COUNTERS;
extern "C" const int const_GL_MAX_GEOMETRY_ATOMIC_COUNTERS;
extern "C" const int const_GL_MAX_FRAGMENT_ATOMIC_COUNTERS;
extern "C" const int const_GL_MAX_COMBINED_ATOMIC_COUNTERS;
extern "C" const int const_GL_MAX_ATOMIC_COUNTER_BUFFER_SIZE;
extern "C" const int const_GL_MAX_ATOMIC_COUNTER_BUFFER_BINDINGS;
extern "C" const int const_GL_ACTIVE_ATOMIC_COUNTER_BUFFERS;
extern "C" const int const_GL_UNIFORM_ATOMIC_COUNTER_BUFFER_INDEX;
extern "C" const int const_GL_UNSIGNED_INT_ATOMIC_COUNTER;
extern "C" const int const_GL_VERTEX_ATTRIB_ARRAY_BARRIER_BIT;
extern "C" const int const_GL_ELEMENT_ARRAY_BARRIER_BIT;
extern "C" const int const_GL_UNIFORM_BARRIER_BIT;
extern "C" const int const_GL_TEXTURE_FETCH_BARRIER_BIT;
extern "C" const int const_GL_SHADER_IMAGE_ACCESS_BARRIER_BIT;
extern "C" const int const_GL_COMMAND_BARRIER_BIT;
extern "C" const int const_GL_PIXEL_BUFFER_BARRIER_BIT;
extern "C" const int const_GL_TEXTURE_UPDATE_BARRIER_BIT;
extern "C" const int const_GL_BUFFER_UPDATE_BARRIER_BIT;
extern "C" const int const_GL_FRAMEBUFFER_BARRIER_BIT;
extern "C" const int const_GL_TRANSFORM_FEEDBACK_BARRIER_BIT;
extern "C" const int const_GL_ATOMIC_COUNTER_BARRIER_BIT;
extern "C" const int const_GL_ALL_BARRIER_BITS;
extern "C" const int const_GL_MAX_IMAGE_UNITS;
extern "C" const int const_GL_MAX_COMBINED_IMAGE_UNITS_AND_FRAGMENT_OUTPUTS;
extern "C" const int const_GL_IMAGE_BINDING_NAME;
extern "C" const int const_GL_IMAGE_BINDING_LEVEL;
extern "C" const int const_GL_IMAGE_BINDING_LAYERED;
extern "C" const int const_GL_IMAGE_BINDING_LAYER;
extern "C" const int const_GL_IMAGE_BINDING_ACCESS;
extern "C" const int const_GL_IMAGE_1D;
extern "C" const int const_GL_IMAGE_2D;
extern "C" const int const_GL_IMAGE_3D;
extern "C" const int const_GL_IMAGE_2D_RECT;
extern "C" const int const_GL_IMAGE_CUBE;
extern "C" const int const_GL_IMAGE_BUFFER;
extern "C" const int const_GL_IMAGE_1D_ARRAY;
extern "C" const int const_GL_IMAGE_2D_ARRAY;
extern "C" const int const_GL_IMAGE_CUBE_MAP_ARRAY;
extern "C" const int const_GL_IMAGE_2D_MULTISAMPLE;
extern "C" const int const_GL_IMAGE_2D_MULTISAMPLE_ARRAY;
extern "C" const int const_GL_INT_IMAGE_1D;
extern "C" const int const_GL_INT_IMAGE_2D;
extern "C" const int const_GL_INT_IMAGE_3D;
extern "C" const int const_GL_INT_IMAGE_2D_RECT;
extern "C" const int const_GL_INT_IMAGE_CUBE;
extern "C" const int const_GL_INT_IMAGE_BUFFER;
extern "C" const int const_GL_INT_IMAGE_1D_ARRAY;
extern "C" const int const_GL_INT_IMAGE_2D_ARRAY;
extern "C" const int const_GL_INT_IMAGE_CUBE_MAP_ARRAY;
extern "C" const int const_GL_INT_IMAGE_2D_MULTISAMPLE;
extern "C" const int const_GL_INT_IMAGE_2D_MULTISAMPLE_ARRAY;
extern "C" const int const_GL_UNSIGNED_INT_IMAGE_1D;
extern "C" const int const_GL_UNSIGNED_INT_IMAGE_2D;
extern "C" const int const_GL_UNSIGNED_INT_IMAGE_3D;
extern "C" const int const_GL_UNSIGNED_INT_IMAGE_2D_RECT;
extern "C" const int const_GL_UNSIGNED_INT_IMAGE_CUBE;
extern "C" const int const_GL_UNSIGNED_INT_IMAGE_BUFFER;
extern "C" const int const_GL_UNSIGNED_INT_IMAGE_1D_ARRAY;
extern "C" const int const_GL_UNSIGNED_INT_IMAGE_2D_ARRAY;
extern "C" const int const_GL_UNSIGNED_INT_IMAGE_CUBE_MAP_ARRAY;
extern "C" const int const_GL_UNSIGNED_INT_IMAGE_2D_MULTISAMPLE;
extern "C" const int const_GL_UNSIGNED_INT_IMAGE_2D_MULTISAMPLE_ARRAY;
extern "C" const int const_GL_MAX_IMAGE_SAMPLES;
extern "C" const int const_GL_IMAGE_BINDING_FORMAT;
extern "C" const int const_GL_IMAGE_FORMAT_COMPATIBILITY_TYPE;
extern "C" const int const_GL_IMAGE_FORMAT_COMPATIBILITY_BY_SIZE;
extern "C" const int const_GL_IMAGE_FORMAT_COMPATIBILITY_BY_CLASS;
extern "C" const int const_GL_MAX_VERTEX_IMAGE_UNIFORMS;
extern "C" const int const_GL_MAX_TESS_CONTROL_IMAGE_UNIFORMS;
extern "C" const int const_GL_MAX_TESS_EVALUATION_IMAGE_UNIFORMS;
extern "C" const int const_GL_MAX_GEOMETRY_IMAGE_UNIFORMS;
extern "C" const int const_GL_MAX_FRAGMENT_IMAGE_UNIFORMS;
extern "C" const int const_GL_MAX_COMBINED_IMAGE_UNIFORMS;
extern "C" const int const_GL_TEXTURE_IMMUTABLE_FORMAT;
extern "C" const int const_GL_COMPRESSED_RGBA_ASTC_4x4_KHR;
extern "C" const int const_GL_COMPRESSED_RGBA_ASTC_5x4_KHR;
extern "C" const int const_GL_COMPRESSED_RGBA_ASTC_5x5_KHR;
extern "C" const int const_GL_COMPRESSED_RGBA_ASTC_6x5_KHR;
extern "C" const int const_GL_COMPRESSED_RGBA_ASTC_6x6_KHR;
extern "C" const int const_GL_COMPRESSED_RGBA_ASTC_8x5_KHR;
extern "C" const int const_GL_COMPRESSED_RGBA_ASTC_8x6_KHR;
extern "C" const int const_GL_COMPRESSED_RGBA_ASTC_8x8_KHR;
extern "C" const int const_GL_COMPRESSED_RGBA_ASTC_10x5_KHR;
extern "C" const int const_GL_COMPRESSED_RGBA_ASTC_10x6_KHR;
extern "C" const int const_GL_COMPRESSED_RGBA_ASTC_10x8_KHR;
extern "C" const int const_GL_COMPRESSED_RGBA_ASTC_10x10_KHR;
extern "C" const int const_GL_COMPRESSED_RGBA_ASTC_12x10_KHR;
extern "C" const int const_GL_COMPRESSED_RGBA_ASTC_12x12_KHR;
extern "C" const int const_GL_COMPRESSED_SRGB8_ALPHA8_ASTC_4x4_KHR;
extern "C" const int const_GL_COMPRESSED_SRGB8_ALPHA8_ASTC_5x4_KHR;
extern "C" const int const_GL_COMPRESSED_SRGB8_ALPHA8_ASTC_5x5_KHR;
extern "C" const int const_GL_COMPRESSED_SRGB8_ALPHA8_ASTC_6x5_KHR;
extern "C" const int const_GL_COMPRESSED_SRGB8_ALPHA8_ASTC_6x6_KHR;
extern "C" const int const_GL_COMPRESSED_SRGB8_ALPHA8_ASTC_8x5_KHR;
extern "C" const int const_GL_COMPRESSED_SRGB8_ALPHA8_ASTC_8x6_KHR;
extern "C" const int const_GL_COMPRESSED_SRGB8_ALPHA8_ASTC_8x8_KHR;
extern "C" const int const_GL_COMPRESSED_SRGB8_ALPHA8_ASTC_10x5_KHR;
extern "C" const int const_GL_COMPRESSED_SRGB8_ALPHA8_ASTC_10x6_KHR;
extern "C" const int const_GL_COMPRESSED_SRGB8_ALPHA8_ASTC_10x8_KHR;
extern "C" const int const_GL_COMPRESSED_SRGB8_ALPHA8_ASTC_10x10_KHR;
extern "C" const int const_GL_COMPRESSED_SRGB8_ALPHA8_ASTC_12x10_KHR;
extern "C" const int const_GL_COMPRESSED_SRGB8_ALPHA8_ASTC_12x12_KHR;
extern "C" const int const_GL_DEBUG_OUTPUT_SYNCHRONOUS;
extern "C" const int const_GL_DEBUG_NEXT_LOGGED_MESSAGE_LENGTH;
extern "C" const int const_GL_DEBUG_CALLBACK_FUNCTION;
extern "C" const int const_GL_DEBUG_CALLBACK_USER_PARAM;
extern "C" const int const_GL_DEBUG_SOURCE_API;
extern "C" const int const_GL_DEBUG_SOURCE_WINDOW_SYSTEM;
extern "C" const int const_GL_DEBUG_SOURCE_SHADER_COMPILER;
extern "C" const int const_GL_DEBUG_SOURCE_THIRD_PARTY;
extern "C" const int const_GL_DEBUG_SOURCE_APPLICATION;
extern "C" const int const_GL_DEBUG_SOURCE_OTHER;
extern "C" const int const_GL_DEBUG_TYPE_ERROR;
extern "C" const int const_GL_DEBUG_TYPE_DEPRECATED_BEHAVIOR;
extern "C" const int const_GL_DEBUG_TYPE_UNDEFINED_BEHAVIOR;
extern "C" const int const_GL_DEBUG_TYPE_PORTABILITY;
extern "C" const int const_GL_DEBUG_TYPE_PERFORMANCE;
extern "C" const int const_GL_DEBUG_TYPE_OTHER;
extern "C" const int const_GL_DEBUG_TYPE_MARKER;
extern "C" const int const_GL_DEBUG_TYPE_PUSH_GROUP;
extern "C" const int const_GL_DEBUG_TYPE_POP_GROUP;
extern "C" const int const_GL_DEBUG_SEVERITY_NOTIFICATION;
extern "C" const int const_GL_MAX_DEBUG_GROUP_STACK_DEPTH;
extern "C" const int const_GL_DEBUG_GROUP_STACK_DEPTH;
extern "C" const int const_GL_BUFFER;
extern "C" const int const_GL_SHADER;
extern "C" const int const_GL_PROGRAM;
extern "C" const int const_GL_QUERY;
extern "C" const int const_GL_PROGRAM_PIPELINE;
extern "C" const int const_GL_SAMPLER;
extern "C" const int const_GL_DISPLAY_LIST;
extern "C" const int const_GL_MAX_LABEL_LENGTH;
extern "C" const int const_GL_MAX_DEBUG_MESSAGE_LENGTH;
extern "C" const int const_GL_MAX_DEBUG_LOGGED_MESSAGES;
extern "C" const int const_GL_DEBUG_LOGGED_MESSAGES;
extern "C" const int const_GL_DEBUG_SEVERITY_HIGH;
extern "C" const int const_GL_DEBUG_SEVERITY_MEDIUM;
extern "C" const int const_GL_DEBUG_SEVERITY_LOW;
extern "C" const int const_GL_DEBUG_OUTPUT;
extern "C" const int const_GL_CONTEXT_FLAG_DEBUG_BIT;
extern "C" const int const_GL_COMPUTE_SHADER;
extern "C" const int const_GL_MAX_COMPUTE_UNIFORM_BLOCKS;
extern "C" const int const_GL_MAX_COMPUTE_TEXTURE_IMAGE_UNITS;
extern "C" const int const_GL_MAX_COMPUTE_IMAGE_UNIFORMS;
extern "C" const int const_GL_MAX_COMPUTE_SHARED_MEMORY_SIZE;
extern "C" const int const_GL_MAX_COMPUTE_UNIFORM_COMPONENTS;
extern "C" const int const_GL_MAX_COMPUTE_ATOMIC_COUNTER_BUFFERS;
extern "C" const int const_GL_MAX_COMPUTE_ATOMIC_COUNTERS;
extern "C" const int const_GL_MAX_COMBINED_COMPUTE_UNIFORM_COMPONENTS;
extern "C" const int const_GL_MAX_COMPUTE_LOCAL_INVOCATIONS;
extern "C" const int const_GL_MAX_COMPUTE_WORK_GROUP_COUNT;
extern "C" const int const_GL_MAX_COMPUTE_WORK_GROUP_SIZE;
extern "C" const int const_GL_COMPUTE_LOCAL_WORK_SIZE;
extern "C" const int const_GL_UNIFORM_BLOCK_REFERENCED_BY_COMPUTE_SHADER;
extern "C" const int const_GL_ATOMIC_COUNTER_BUFFER_REFERENCED_BY_COMPUTE_SHADER;
extern "C" const int const_GL_DISPATCH_INDIRECT_BUFFER;
extern "C" const int const_GL_DISPATCH_INDIRECT_BUFFER_BINDING;
extern "C" const int const_GL_COMPUTE_SHADER_BIT;
extern "C" const int const_GL_TEXTURE_VIEW_MIN_LEVEL;
extern "C" const int const_GL_TEXTURE_VIEW_NUM_LEVELS;
extern "C" const int const_GL_TEXTURE_VIEW_MIN_LAYER;
extern "C" const int const_GL_TEXTURE_VIEW_NUM_LAYERS;
extern "C" const int const_GL_TEXTURE_IMMUTABLE_LEVELS;
extern "C" const int const_GL_VERTEX_ATTRIB_BINDING;
extern "C" const int const_GL_VERTEX_ATTRIB_RELATIVE_OFFSET;
extern "C" const int const_GL_VERTEX_BINDING_DIVISOR;
extern "C" const int const_GL_VERTEX_BINDING_OFFSET;
extern "C" const int const_GL_VERTEX_BINDING_STRIDE;
extern "C" const int const_GL_MAX_VERTEX_ATTRIB_RELATIVE_OFFSET;
extern "C" const int const_GL_MAX_VERTEX_ATTRIB_BINDINGS;
extern "C" const int const_GL_COMPRESSED_RGB8_ETC2;
extern "C" const int const_GL_COMPRESSED_SRGB8_ETC2;
extern "C" const int const_GL_COMPRESSED_RGB8_PUNCHTHROUGH_ALPHA1_ETC2;
extern "C" const int const_GL_COMPRESSED_SRGB8_PUNCHTHROUGH_ALPHA1_ETC2;
extern "C" const int const_GL_COMPRESSED_RGBA8_ETC2_EAC;
extern "C" const int const_GL_COMPRESSED_SRGB8_ALPHA8_ETC2_EAC;
extern "C" const int const_GL_COMPRESSED_R11_EAC;
extern "C" const int const_GL_COMPRESSED_SIGNED_R11_EAC;
extern "C" const int const_GL_COMPRESSED_RG11_EAC;
extern "C" const int const_GL_COMPRESSED_SIGNED_RG11_EAC;
extern "C" const int const_GL_PRIMITIVE_RESTART_FIXED_INDEX;
extern "C" const int const_GL_ANY_SAMPLES_PASSED_CONSERVATIVE;
extern "C" const int const_GL_MAX_ELEMENT_INDEX;
extern "C" const int const_GL_MAX_UNIFORM_LOCATIONS;
extern "C" const int const_GL_FRAMEBUFFER_DEFAULT_WIDTH;
extern "C" const int const_GL_FRAMEBUFFER_DEFAULT_HEIGHT;
extern "C" const int const_GL_FRAMEBUFFER_DEFAULT_LAYERS;
extern "C" const int const_GL_FRAMEBUFFER_DEFAULT_SAMPLES;
extern "C" const int const_GL_FRAMEBUFFER_DEFAULT_FIXED_SAMPLE_LOCATIONS;
extern "C" const int const_GL_MAX_FRAMEBUFFER_WIDTH;
extern "C" const int const_GL_MAX_FRAMEBUFFER_HEIGHT;
extern "C" const int const_GL_MAX_FRAMEBUFFER_LAYERS;
extern "C" const int const_GL_MAX_FRAMEBUFFER_SAMPLES;
extern "C" const int const_GL_INTERNALFORMAT_SUPPORTED;
extern "C" const int const_GL_INTERNALFORMAT_PREFERRED;
extern "C" const int const_GL_INTERNALFORMAT_RED_SIZE;
extern "C" const int const_GL_INTERNALFORMAT_GREEN_SIZE;
extern "C" const int const_GL_INTERNALFORMAT_BLUE_SIZE;
extern "C" const int const_GL_INTERNALFORMAT_ALPHA_SIZE;
extern "C" const int const_GL_INTERNALFORMAT_DEPTH_SIZE;
extern "C" const int const_GL_INTERNALFORMAT_STENCIL_SIZE;
extern "C" const int const_GL_INTERNALFORMAT_SHARED_SIZE;
extern "C" const int const_GL_INTERNALFORMAT_RED_TYPE;
extern "C" const int const_GL_INTERNALFORMAT_GREEN_TYPE;
extern "C" const int const_GL_INTERNALFORMAT_BLUE_TYPE;
extern "C" const int const_GL_INTERNALFORMAT_ALPHA_TYPE;
extern "C" const int const_GL_INTERNALFORMAT_DEPTH_TYPE;
extern "C" const int const_GL_INTERNALFORMAT_STENCIL_TYPE;
extern "C" const int const_GL_MAX_WIDTH;
extern "C" const int const_GL_MAX_HEIGHT;
extern "C" const int const_GL_MAX_DEPTH;
extern "C" const int const_GL_MAX_LAYERS;
extern "C" const int const_GL_MAX_COMBINED_DIMENSIONS;
extern "C" const int const_GL_COLOR_COMPONENTS;
extern "C" const int const_GL_DEPTH_COMPONENTS;
extern "C" const int const_GL_STENCIL_COMPONENTS;
extern "C" const int const_GL_COLOR_RENDERABLE;
extern "C" const int const_GL_DEPTH_RENDERABLE;
extern "C" const int const_GL_STENCIL_RENDERABLE;
extern "C" const int const_GL_FRAMEBUFFER_RENDERABLE;
extern "C" const int const_GL_FRAMEBUFFER_RENDERABLE_LAYERED;
extern "C" const int const_GL_FRAMEBUFFER_BLEND;
extern "C" const int const_GL_READ_PIXELS;
extern "C" const int const_GL_READ_PIXELS_FORMAT;
extern "C" const int const_GL_READ_PIXELS_TYPE;
extern "C" const int const_GL_TEXTURE_IMAGE_FORMAT;
extern "C" const int const_GL_TEXTURE_IMAGE_TYPE;
extern "C" const int const_GL_GET_TEXTURE_IMAGE_FORMAT;
extern "C" const int const_GL_GET_TEXTURE_IMAGE_TYPE;
extern "C" const int const_GL_MIPMAP;
extern "C" const int const_GL_MANUAL_GENERATE_MIPMAP;
extern "C" const int const_GL_AUTO_GENERATE_MIPMAP;
extern "C" const int const_GL_COLOR_ENCODING;
extern "C" const int const_GL_SRGB_READ;
extern "C" const int const_GL_SRGB_WRITE;
extern "C" const int const_GL_SRGB_DECODE_ARB;
extern "C" const int const_GL_FILTER;
extern "C" const int const_GL_VERTEX_TEXTURE;
extern "C" const int const_GL_TESS_CONTROL_TEXTURE;
extern "C" const int const_GL_TESS_EVALUATION_TEXTURE;
extern "C" const int const_GL_GEOMETRY_TEXTURE;
extern "C" const int const_GL_FRAGMENT_TEXTURE;
extern "C" const int const_GL_COMPUTE_TEXTURE;
extern "C" const int const_GL_TEXTURE_SHADOW;
extern "C" const int const_GL_TEXTURE_GATHER;
extern "C" const int const_GL_TEXTURE_GATHER_SHADOW;
extern "C" const int const_GL_SHADER_IMAGE_LOAD;
extern "C" const int const_GL_SHADER_IMAGE_STORE;
extern "C" const int const_GL_SHADER_IMAGE_ATOMIC;
extern "C" const int const_GL_IMAGE_TEXEL_SIZE;
extern "C" const int const_GL_IMAGE_COMPATIBILITY_CLASS;
extern "C" const int const_GL_IMAGE_PIXEL_FORMAT;
extern "C" const int const_GL_IMAGE_PIXEL_TYPE;
extern "C" const int const_GL_SIMULTANEOUS_TEXTURE_AND_DEPTH_TEST;
extern "C" const int const_GL_SIMULTANEOUS_TEXTURE_AND_STENCIL_TEST;
extern "C" const int const_GL_SIMULTANEOUS_TEXTURE_AND_DEPTH_WRITE;
extern "C" const int const_GL_SIMULTANEOUS_TEXTURE_AND_STENCIL_WRITE;
extern "C" const int const_GL_TEXTURE_COMPRESSED_BLOCK_WIDTH;
extern "C" const int const_GL_TEXTURE_COMPRESSED_BLOCK_HEIGHT;
extern "C" const int const_GL_TEXTURE_COMPRESSED_BLOCK_SIZE;
extern "C" const int const_GL_CLEAR_BUFFER;
extern "C" const int const_GL_TEXTURE_VIEW;
extern "C" const int const_GL_VIEW_COMPATIBILITY_CLASS;
extern "C" const int const_GL_FULL_SUPPORT;
extern "C" const int const_GL_CAVEAT_SUPPORT;
extern "C" const int const_GL_IMAGE_CLASS_4_X_32;
extern "C" const int const_GL_IMAGE_CLASS_2_X_32;
extern "C" const int const_GL_IMAGE_CLASS_1_X_32;
extern "C" const int const_GL_IMAGE_CLASS_4_X_16;
extern "C" const int const_GL_IMAGE_CLASS_2_X_16;
extern "C" const int const_GL_IMAGE_CLASS_1_X_16;
extern "C" const int const_GL_IMAGE_CLASS_4_X_8;
extern "C" const int const_GL_IMAGE_CLASS_2_X_8;
extern "C" const int const_GL_IMAGE_CLASS_1_X_8;
extern "C" const int const_GL_IMAGE_CLASS_11_11_10;
extern "C" const int const_GL_IMAGE_CLASS_10_10_10_2;
extern "C" const int const_GL_VIEW_CLASS_128_BITS;
extern "C" const int const_GL_VIEW_CLASS_96_BITS;
extern "C" const int const_GL_VIEW_CLASS_64_BITS;
extern "C" const int const_GL_VIEW_CLASS_48_BITS;
extern "C" const int const_GL_VIEW_CLASS_32_BITS;
extern "C" const int const_GL_VIEW_CLASS_24_BITS;
extern "C" const int const_GL_VIEW_CLASS_16_BITS;
extern "C" const int const_GL_VIEW_CLASS_8_BITS;
extern "C" const int const_GL_VIEW_CLASS_S3TC_DXT1_RGB;
extern "C" const int const_GL_VIEW_CLASS_S3TC_DXT1_RGBA;
extern "C" const int const_GL_VIEW_CLASS_S3TC_DXT3_RGBA;
extern "C" const int const_GL_VIEW_CLASS_S3TC_DXT5_RGBA;
extern "C" const int const_GL_VIEW_CLASS_RGTC1_RED;
extern "C" const int const_GL_VIEW_CLASS_RGTC2_RG;
extern "C" const int const_GL_VIEW_CLASS_BPTC_UNORM;
extern "C" const int const_GL_VIEW_CLASS_BPTC_FLOAT;
extern "C" const int const_GL_UNIFORM;
extern "C" const int const_GL_UNIFORM_BLOCK;
extern "C" const int const_GL_PROGRAM_INPUT;
extern "C" const int const_GL_PROGRAM_OUTPUT;
extern "C" const int const_GL_BUFFER_VARIABLE;
extern "C" const int const_GL_SHADER_STORAGE_BLOCK;
extern "C" const int const_GL_VERTEX_SUBROUTINE;
extern "C" const int const_GL_TESS_CONTROL_SUBROUTINE;
extern "C" const int const_GL_TESS_EVALUATION_SUBROUTINE;
extern "C" const int const_GL_GEOMETRY_SUBROUTINE;
extern "C" const int const_GL_FRAGMENT_SUBROUTINE;
extern "C" const int const_GL_COMPUTE_SUBROUTINE;
extern "C" const int const_GL_VERTEX_SUBROUTINE_UNIFORM;
extern "C" const int const_GL_TESS_CONTROL_SUBROUTINE_UNIFORM;
extern "C" const int const_GL_TESS_EVALUATION_SUBROUTINE_UNIFORM;
extern "C" const int const_GL_GEOMETRY_SUBROUTINE_UNIFORM;
extern "C" const int const_GL_FRAGMENT_SUBROUTINE_UNIFORM;
extern "C" const int const_GL_COMPUTE_SUBROUTINE_UNIFORM;
extern "C" const int const_GL_TRANSFORM_FEEDBACK_VARYING;
extern "C" const int const_GL_ACTIVE_RESOURCES;
extern "C" const int const_GL_MAX_NAME_LENGTH;
extern "C" const int const_GL_MAX_NUM_ACTIVE_VARIABLES;
extern "C" const int const_GL_MAX_NUM_COMPATIBLE_SUBROUTINES;
extern "C" const int const_GL_NAME_LENGTH;
extern "C" const int const_GL_TYPE;
extern "C" const int const_GL_ARRAY_SIZE;
extern "C" const int const_GL_OFFSET;
extern "C" const int const_GL_BLOCK_INDEX;
extern "C" const int const_GL_ARRAY_STRIDE;
extern "C" const int const_GL_MATRIX_STRIDE;
extern "C" const int const_GL_IS_ROW_MAJOR;
extern "C" const int const_GL_ATOMIC_COUNTER_BUFFER_INDEX;
extern "C" const int const_GL_BUFFER_BINDING;
extern "C" const int const_GL_BUFFER_DATA_SIZE;
extern "C" const int const_GL_NUM_ACTIVE_VARIABLES;
extern "C" const int const_GL_ACTIVE_VARIABLES;
extern "C" const int const_GL_REFERENCED_BY_VERTEX_SHADER;
extern "C" const int const_GL_REFERENCED_BY_TESS_CONTROL_SHADER;
extern "C" const int const_GL_REFERENCED_BY_TESS_EVALUATION_SHADER;
extern "C" const int const_GL_REFERENCED_BY_GEOMETRY_SHADER;
extern "C" const int const_GL_REFERENCED_BY_FRAGMENT_SHADER;
extern "C" const int const_GL_REFERENCED_BY_COMPUTE_SHADER;
extern "C" const int const_GL_TOP_LEVEL_ARRAY_SIZE;
extern "C" const int const_GL_TOP_LEVEL_ARRAY_STRIDE;
extern "C" const int const_GL_LOCATION;
extern "C" const int const_GL_LOCATION_INDEX;
extern "C" const int const_GL_IS_PER_PATCH;
extern "C" const int const_GL_SHADER_STORAGE_BUFFER;
extern "C" const int const_GL_SHADER_STORAGE_BUFFER_BINDING;
extern "C" const int const_GL_SHADER_STORAGE_BUFFER_START;
extern "C" const int const_GL_SHADER_STORAGE_BUFFER_SIZE;
extern "C" const int const_GL_MAX_VERTEX_SHADER_STORAGE_BLOCKS;
extern "C" const int const_GL_MAX_GEOMETRY_SHADER_STORAGE_BLOCKS;
extern "C" const int const_GL_MAX_TESS_CONTROL_SHADER_STORAGE_BLOCKS;
extern "C" const int const_GL_MAX_TESS_EVALUATION_SHADER_STORAGE_BLOCKS;
extern "C" const int const_GL_MAX_FRAGMENT_SHADER_STORAGE_BLOCKS;
extern "C" const int const_GL_MAX_COMPUTE_SHADER_STORAGE_BLOCKS;
extern "C" const int const_GL_MAX_COMBINED_SHADER_STORAGE_BLOCKS;
extern "C" const int const_GL_MAX_SHADER_STORAGE_BUFFER_BINDINGS;
extern "C" const int const_GL_MAX_SHADER_STORAGE_BLOCK_SIZE;
extern "C" const int const_GL_SHADER_STORAGE_BUFFER_OFFSET_ALIGNMENT;
extern "C" const int const_GL_SHADER_STORAGE_BARRIER_BIT;
extern "C" const int const_GL_MAX_COMBINED_SHADER_OUTPUT_RESOURCES;
extern "C" const int const_GL_DEPTH_STENCIL_TEXTURE_MODE;
extern "C" const int const_GL_TEXTURE_BUFFER_OFFSET;
extern "C" const int const_GL_TEXTURE_BUFFER_SIZE;
extern "C" const int const_GL_TEXTURE_BUFFER_OFFSET_ALIGNMENT;
extern "C" const int const_GLEXT_64_TYPES_DEFINED;
extern "C" const int const_GL_VERSION_1_0;
extern "C" const int const_GL_VERSION_1_1;
extern "C" const int const_GL_VERSION_1_2;
extern "C" const int const_GL_VERSION_1_3;
extern "C" const int const_GL_VERSION_1_4;
extern "C" const int const_GL_VERSION_1_5;
extern "C" const int const_GL_VERSION_2_0;
extern "C" const int const_GL_VERSION_2_1;
extern "C" const int const_GL_VERSION_3_0;
extern "C" const int const_GL_VERSION_3_1;
extern "C" const int const_GL_VERSION_3_2;
extern "C" const int const_GL_VERSION_3_3;
extern "C" const int const_GL_VERSION_4_0;
extern "C" const int const_GL_VERSION_4_1;
extern "C" const int const_GL_VERSION_4_2;
extern "C" const int const_GL_VERSION_4_3;
extern "C" const int const_GL_ARB_depth_buffer_float;
extern "C" const int const_GL_ARB_framebuffer_object;
extern "C" const int const_GL_ARB_framebuffer_sRGB;
extern "C" const int const_GL_ARB_half_float_vertex;
extern "C" const int const_GL_ARB_map_buffer_range;
extern "C" const int const_GL_ARB_texture_compression_rgtc;
extern "C" const int const_GL_ARB_texture_rg;
extern "C" const int const_GL_ARB_vertex_array_object;
extern "C" const int const_GL_ARB_uniform_buffer_object;
extern "C" const int const_GL_ARB_copy_buffer;
extern "C" const int const_GL_ARB_depth_clamp;
extern "C" const int const_GL_ARB_draw_elements_base_vertex;
extern "C" const int const_GL_ARB_fragment_coord_conventions;
extern "C" const int const_GL_ARB_provoking_vertex;
extern "C" const int const_GL_ARB_seamless_cube_map;
extern "C" const int const_GL_ARB_sync;
extern "C" const int const_GL_ARB_texture_multisample;
extern "C" const int const_GL_ARB_vertex_array_bgra;
extern "C" const int const_GL_ARB_draw_buffers_blend;
extern "C" const int const_GL_ARB_sample_shading;
extern "C" const int const_GL_ARB_texture_cube_map_array;
extern "C" const int const_GL_ARB_texture_gather;
extern "C" const int const_GL_ARB_texture_query_lod;
extern "C" const int const_GL_ARB_shading_language_include;
extern "C" const int const_GL_ARB_texture_compression_bptc;
extern "C" const int const_GL_ARB_blend_func_extended;
extern "C" const int const_GL_ARB_explicit_attrib_location;
extern "C" const int const_GL_ARB_occlusion_query2;
extern "C" const int const_GL_ARB_sampler_objects;
extern "C" const int const_GL_ARB_shader_bit_encoding;
extern "C" const int const_GL_ARB_texture_rgb10_a2ui;
extern "C" const int const_GL_ARB_texture_swizzle;
extern "C" const int const_GL_ARB_timer_query;
extern "C" const int const_GL_ARB_vertex_type_2_10_10_10_rev;
extern "C" const int const_GL_ARB_draw_indirect;
extern "C" const int const_GL_ARB_gpu_shader5;
extern "C" const int const_GL_ARB_gpu_shader_fp64;
extern "C" const int const_GL_ARB_shader_subroutine;
extern "C" const int const_GL_ARB_tessellation_shader;
extern "C" const int const_GL_ARB_texture_buffer_object_rgb32;
extern "C" const int const_GL_ARB_transform_feedback2;
extern "C" const int const_GL_ARB_transform_feedback3;
extern "C" const int const_GL_ARB_ES2_compatibility;
extern "C" const int const_GL_ARB_get_program_binary;
extern "C" const int const_GL_ARB_separate_shader_objects;
extern "C" const int const_GL_ARB_vertex_attrib_64bit;
extern "C" const int const_GL_ARB_viewport_array;
extern "C" const int const_GL_ARB_cl_event;
extern "C" const int const_GL_ARB_debug_output;
extern "C" const int const_GL_ARB_robustness;
extern "C" const int const_GL_ARB_shader_stencil_export;
extern "C" const int const_GL_ARB_base_instance;
extern "C" const int const_GL_ARB_shading_language_420pack;
extern "C" const int const_GL_ARB_transform_feedback_instanced;
extern "C" const int const_GL_ARB_compressed_texture_pixel_storage;
extern "C" const int const_GL_ARB_conservative_depth;
extern "C" const int const_GL_ARB_internalformat_query;
extern "C" const int const_GL_ARB_map_buffer_alignment;
extern "C" const int const_GL_ARB_shader_atomic_counters;
extern "C" const int const_GL_ARB_shader_image_load_store;
extern "C" const int const_GL_ARB_shading_language_packing;
extern "C" const int const_GL_ARB_texture_storage;
extern "C" const int const_GL_KHR_texture_compression_astc_ldr;
extern "C" const int const_GL_KHR_debug;
extern "C" const int const_GL_ARB_arrays_of_arrays;
extern "C" const int const_GL_ARB_clear_buffer_object;
extern "C" const int const_GL_ARB_compute_shader;
extern "C" const int const_GL_ARB_copy_image;
extern "C" const int const_GL_ARB_texture_view;
extern "C" const int const_GL_ARB_vertex_attrib_binding;
extern "C" const int const_GL_ARB_robustness_isolation;
extern "C" const int const_GL_ARB_ES3_compatibility;
extern "C" const int const_GL_ARB_explicit_uniform_location;
extern "C" const int const_GL_ARB_fragment_layer_viewport;
extern "C" const int const_GL_ARB_framebuffer_no_attachments;
extern "C" const int const_GL_ARB_internalformat_query2;
extern "C" const int const_GL_ARB_invalidate_subdata;
extern "C" const int const_GL_ARB_multi_draw_indirect;
extern "C" const int const_GL_ARB_program_interface_query;
extern "C" const int const_GL_ARB_robust_buffer_access_behavior;
extern "C" const int const_GL_ARB_shader_image_size;
extern "C" const int const_GL_ARB_shader_storage_buffer_object;
extern "C" const int const_GL_ARB_stencil_texturing;
extern "C" const int const_GL_ARB_texture_buffer_range;
extern "C" const int const_GL_ARB_texture_query_levels;
extern "C" const int const_GL_ARB_texture_storage_multisample;

extern "C" // GLFW;
extern "C" const int const_GLFW_TRUE;
extern "C" const int const_GLFW_FALSE;
extern "C" const int const_GLFW_SAMPLES;
extern "C" const int const_GLFW_CONTEXT_VERSION_MAJOR;
extern "C" const int const_GLFW_CONTEXT_VERSION_MINOR;
extern "C" const int const_GLFW_OPENGL_PROFILE;
extern "C" const int const_GLFW_OPENGL_FORWARD_COMPAT;
extern "C" const int const_GLFW_OPENGL_CORE_PROFILE;
extern "C" const int const_GLFW_STICKY_KEYS;

// cimgui
// enum ImGuiWindowFlags_
extern "C" const int const_ImGuiWindowFlags_None;
extern "C" const int const_ImGuiWindowFlags_NoTitleBar;
extern "C" const int const_ImGuiWindowFlags_NoResize;
extern "C" const int const_ImGuiWindowFlags_NoMove;
extern "C" const int const_ImGuiWindowFlags_NoScrollbar;
extern "C" const int const_ImGuiWindowFlags_NoScrollWithMouse;
extern "C" const int const_ImGuiWindowFlags_NoCollapse;
extern "C" const int const_ImGuiWindowFlags_AlwaysAutoResize;
extern "C" const int const_ImGuiWindowFlags_NoBackground;
extern "C" const int const_ImGuiWindowFlags_NoSavedSettings;
extern "C" const int const_ImGuiWindowFlags_NoMouseInputs;
extern "C" const int const_ImGuiWindowFlags_MenuBar;
extern "C" const int const_ImGuiWindowFlags_HorizontalScrollbar;
extern "C" const int const_ImGuiWindowFlags_NoFocusOnAppearing;
extern "C" const int const_ImGuiWindowFlags_NoBringToFrontOnFocus;
extern "C" const int const_ImGuiWindowFlags_AlwaysVerticalScrollbar;
extern "C" const int const_ImGuiWindowFlags_AlwaysHorizontalScrollbar;
extern "C" const int const_ImGuiWindowFlags_AlwaysUseWindowPadding;
extern "C" const int const_ImGuiWindowFlags_NoNavInputs;
extern "C" const int const_ImGuiWindowFlags_NoNavFocus;
extern "C" const int const_ImGuiWindowFlags_UnsavedDocument;
extern "C" const int const_ImGuiWindowFlags_NoDocking;
extern "C" const int const_ImGuiWindowFlags_NoNav;
extern "C" const int const_ImGuiWindowFlags_NoDecoration;
extern "C" const int const_ImGuiWindowFlags_NoInputs;
extern "C" const int const_ImGuiWindowFlags_NavFlattened;
extern "C" const int const_ImGuiWindowFlags_ChildWindow;
extern "C" const int const_ImGuiWindowFlags_Tooltip;
extern "C" const int const_ImGuiWindowFlags_Popup;
extern "C" const int const_ImGuiWindowFlags_Modal;
extern "C" const int const_ImGuiWindowFlags_ChildMenu;
extern "C" const int const_ImGuiWindowFlags_DockNodeHost;
// enum ImGuiInputTextFlags_
extern "C" const int const_ImGuiInputTextFlags_None;
extern "C" const int const_ImGuiInputTextFlags_CharsDecimal;
extern "C" const int const_ImGuiInputTextFlags_CharsHexadecimal;
extern "C" const int const_ImGuiInputTextFlags_CharsUppercase;
extern "C" const int const_ImGuiInputTextFlags_CharsNoBlank;
extern "C" const int const_ImGuiInputTextFlags_AutoSelectAll;
extern "C" const int const_ImGuiInputTextFlags_EnterReturnsTrue;
extern "C" const int const_ImGuiInputTextFlags_CallbackCompletion;
extern "C" const int const_ImGuiInputTextFlags_CallbackHistory;
extern "C" const int const_ImGuiInputTextFlags_CallbackAlways;
extern "C" const int const_ImGuiInputTextFlags_CallbackCharFilter;
extern "C" const int const_ImGuiInputTextFlags_AllowTabInput;
extern "C" const int const_ImGuiInputTextFlags_CtrlEnterForNewLine;
extern "C" const int const_ImGuiInputTextFlags_NoHorizontalScroll;
extern "C" const int const_ImGuiInputTextFlags_AlwaysOverwrite;
extern "C" const int const_ImGuiInputTextFlags_ReadOnly;
extern "C" const int const_ImGuiInputTextFlags_Password;
extern "C" const int const_ImGuiInputTextFlags_NoUndoRedo;
extern "C" const int const_ImGuiInputTextFlags_CharsScientific;
extern "C" const int const_ImGuiInputTextFlags_CallbackResize;
extern "C" const int const_ImGuiInputTextFlags_CallbackEdit;
// enum ImGuiTreeNodeFlags_;
extern "C" const int const_ImGuiTreeNodeFlags_None;
extern "C" const int const_ImGuiTreeNodeFlags_Selected;
extern "C" const int const_ImGuiTreeNodeFlags_Framed;
extern "C" const int const_ImGuiTreeNodeFlags_AllowItemOverlap;
extern "C" const int const_ImGuiTreeNodeFlags_NoTreePushOnOpen;
extern "C" const int const_ImGuiTreeNodeFlags_NoAutoOpenOnLog;
extern "C" const int const_ImGuiTreeNodeFlags_DefaultOpen;
extern "C" const int const_ImGuiTreeNodeFlags_OpenOnDoubleClick;
extern "C" const int const_ImGuiTreeNodeFlags_OpenOnArrow;
extern "C" const int const_ImGuiTreeNodeFlags_Leaf;
extern "C" const int const_ImGuiTreeNodeFlags_Bullet;
extern "C" const int const_ImGuiTreeNodeFlags_FramePadding;
extern "C" const int const_ImGuiTreeNodeFlags_SpanAvailWidth;
extern "C" const int const_ImGuiTreeNodeFlags_SpanFullWidth;
extern "C" const int const_ImGuiTreeNodeFlags_NavLeftJumpsBackHere;
extern "C" const int const_ImGuiTreeNodeFlags_CollapsingHeader;
// enum ImGuiPopupFlags_
extern "C" const int const_ImGuiPopupFlags_None;
extern "C" const int const_ImGuiPopupFlags_MouseButtonLeft;
extern "C" const int const_ImGuiPopupFlags_MouseButtonRight;
extern "C" const int const_ImGuiPopupFlags_MouseButtonMiddle;
extern "C" const int const_ImGuiPopupFlags_MouseButtonMask_;
extern "C" const int const_ImGuiPopupFlags_MouseButtonDefault_;
extern "C" const int const_ImGuiPopupFlags_NoOpenOverExistingPopup;
extern "C" const int const_ImGuiPopupFlags_NoOpenOverItems;
extern "C" const int const_ImGuiPopupFlags_AnyPopupId;
extern "C" const int const_ImGuiPopupFlags_AnyPopupLevel;
extern "C" const int const_ImGuiPopupFlags_AnyPopup;
// enum ImGuiSelectableFlags_
extern "C" const int const_ImGuiSelectableFlags_None;
extern "C" const int const_ImGuiSelectableFlags_DontClosePopups;
extern "C" const int const_ImGuiSelectableFlags_SpanAllColumns;
extern "C" const int const_ImGuiSelectableFlags_AllowDoubleClick;
extern "C" const int const_ImGuiSelectableFlags_Disabled;
extern "C" const int const_ImGuiSelectableFlags_AllowItemOverlap;
// enum ImGuiComboFlags_
extern "C" const int const_ImGuiComboFlags_None;
extern "C" const int const_ImGuiComboFlags_PopupAlignLeft;
extern "C" const int const_ImGuiComboFlags_HeightSmall;
extern "C" const int const_ImGuiComboFlags_HeightRegular;
extern "C" const int const_ImGuiComboFlags_HeightLarge;
extern "C" const int const_ImGuiComboFlags_HeightLargest;
extern "C" const int const_ImGuiComboFlags_NoArrowButton;
extern "C" const int const_ImGuiComboFlags_NoPreview;
extern "C" const int const_ImGuiComboFlags_HeightMask_;
// enum ImGuiTabBarFlags_
extern "C" const int const_ImGuiTabBarFlags_None;
extern "C" const int const_ImGuiTabBarFlags_Reorderable;
extern "C" const int const_ImGuiTabBarFlags_AutoSelectNewTabs;
extern "C" const int const_ImGuiTabBarFlags_TabListPopupButton;
extern "C" const int const_ImGuiTabBarFlags_NoCloseWithMiddleMouseButton;
extern "C" const int const_ImGuiTabBarFlags_NoTabListScrollingButtons;
extern "C" const int const_ImGuiTabBarFlags_NoTooltip;
extern "C" const int const_ImGuiTabBarFlags_FittingPolicyResizeDown;
extern "C" const int const_ImGuiTabBarFlags_FittingPolicyScroll;
extern "C" const int const_ImGuiTabBarFlags_FittingPolicyMask_;
extern "C" const int const_ImGuiTabBarFlags_FittingPolicyDefault_;
// enum ImGuiTabItemFlags_
extern "C" const int const_ImGuiTabItemFlags_None;
extern "C" const int const_ImGuiTabItemFlags_UnsavedDocument;
extern "C" const int const_ImGuiTabItemFlags_SetSelected;
extern "C" const int const_ImGuiTabItemFlags_NoCloseWithMiddleMouseButton;
extern "C" const int const_ImGuiTabItemFlags_NoPushId;
extern "C" const int const_ImGuiTabItemFlags_NoTooltip;
extern "C" const int const_ImGuiTabItemFlags_NoReorder;
extern "C" const int const_ImGuiTabItemFlags_Leading;
extern "C" const int const_ImGuiTabItemFlags_Trailing;
// enum ImGuiTableFlags_
extern "C" const int const_ImGuiTableFlags_None;
extern "C" const int const_ImGuiTableFlags_Resizable;
extern "C" const int const_ImGuiTableFlags_Reorderable;
extern "C" const int const_ImGuiTableFlags_Hideable;
extern "C" const int const_ImGuiTableFlags_Sortable;
extern "C" const int const_ImGuiTableFlags_NoSavedSettings;
extern "C" const int const_ImGuiTableFlags_ContextMenuInBody;
extern "C" const int const_ImGuiTableFlags_RowBg;
extern "C" const int const_ImGuiTableFlags_BordersInnerH;
extern "C" const int const_ImGuiTableFlags_BordersOuterH;
extern "C" const int const_ImGuiTableFlags_BordersInnerV;
extern "C" const int const_ImGuiTableFlags_BordersOuterV;
extern "C" const int const_ImGuiTableFlags_BordersH;
extern "C" const int const_ImGuiTableFlags_BordersV;
extern "C" const int const_ImGuiTableFlags_BordersInner;
extern "C" const int const_ImGuiTableFlags_BordersOuter;
extern "C" const int const_ImGuiTableFlags_Borders;
extern "C" const int const_ImGuiTableFlags_NoBordersInBody;
extern "C" const int const_ImGuiTableFlags_NoBordersInBodyUntilResize;
extern "C" const int const_ImGuiTableFlags_SizingFixedFit;
extern "C" const int const_ImGuiTableFlags_SizingFixedSame;
extern "C" const int const_ImGuiTableFlags_SizingStretchProp;
extern "C" const int const_ImGuiTableFlags_SizingStretchSame;
extern "C" const int const_ImGuiTableFlags_NoHostExtendX;
extern "C" const int const_ImGuiTableFlags_NoHostExtendY;
extern "C" const int const_ImGuiTableFlags_NoKeepColumnsVisible;
extern "C" const int const_ImGuiTableFlags_PreciseWidths;
extern "C" const int const_ImGuiTableFlags_NoClip;
extern "C" const int const_ImGuiTableFlags_PadOuterX;
extern "C" const int const_ImGuiTableFlags_NoPadOuterX;
extern "C" const int const_ImGuiTableFlags_NoPadInnerX;
extern "C" const int const_ImGuiTableFlags_ScrollX;
extern "C" const int const_ImGuiTableFlags_ScrollY;
extern "C" const int const_ImGuiTableFlags_SortMulti;
extern "C" const int const_ImGuiTableFlags_SortTristate;
extern "C" const int const_ImGuiTableFlags_SizingMask_;
// enum ImGuiTableColumnFlags_
extern "C" const int const_ImGuiTableColumnFlags_None;
extern "C" const int const_ImGuiTableColumnFlags_Disabled;
extern "C" const int const_ImGuiTableColumnFlags_DefaultHide;
extern "C" const int const_ImGuiTableColumnFlags_DefaultSort;
extern "C" const int const_ImGuiTableColumnFlags_WidthStretch;
extern "C" const int const_ImGuiTableColumnFlags_WidthFixed;
extern "C" const int const_ImGuiTableColumnFlags_NoResize;
extern "C" const int const_ImGuiTableColumnFlags_NoReorder;
extern "C" const int const_ImGuiTableColumnFlags_NoHide;
extern "C" const int const_ImGuiTableColumnFlags_NoClip;
extern "C" const int const_ImGuiTableColumnFlags_NoSort;
extern "C" const int const_ImGuiTableColumnFlags_NoSortAscending;
extern "C" const int const_ImGuiTableColumnFlags_NoSortDescending;
extern "C" const int const_ImGuiTableColumnFlags_NoHeaderLabel;
extern "C" const int const_ImGuiTableColumnFlags_NoHeaderWidth;
extern "C" const int const_ImGuiTableColumnFlags_PreferSortAscending;
extern "C" const int const_ImGuiTableColumnFlags_PreferSortDescending;
extern "C" const int const_ImGuiTableColumnFlags_IndentEnable;
extern "C" const int const_ImGuiTableColumnFlags_IndentDisable;
extern "C" const int const_ImGuiTableColumnFlags_IsEnabled;
extern "C" const int const_ImGuiTableColumnFlags_IsVisible;
extern "C" const int const_ImGuiTableColumnFlags_IsSorted;
extern "C" const int const_ImGuiTableColumnFlags_IsHovered;
extern "C" const int const_ImGuiTableColumnFlags_WidthMask_;
extern "C" const int const_ImGuiTableColumnFlags_IndentMask_;
extern "C" const int const_ImGuiTableColumnFlags_StatusMask_;
extern "C" const int const_ImGuiTableColumnFlags_NoDirectResize_;
// enum ImGuiTableRowFlags_
extern "C" const int const_ImGuiTableRowFlags_None;
extern "C" const int const_ImGuiTableRowFlags_Headers;
// enum ImGuiTableBgTarget_
extern "C" const int const_ImGuiTableBgTarget_None;
extern "C" const int const_ImGuiTableBgTarget_RowBg0;
extern "C" const int const_ImGuiTableBgTarget_RowBg1;
extern "C" const int const_ImGuiTableBgTarget_CellBg;
// enum ImGuiFocusedFlags_
extern "C" const int const_ImGuiFocusedFlags_None;
extern "C" const int const_ImGuiFocusedFlags_ChildWindows;
extern "C" const int const_ImGuiFocusedFlags_RootWindow;
extern "C" const int const_ImGuiFocusedFlags_AnyWindow;
extern "C" const int const_ImGuiFocusedFlags_NoPopupHierarchy;
extern "C" const int const_ImGuiFocusedFlags_DockHierarchy;
extern "C" const int const_ImGuiFocusedFlags_RootAndChildWindows;
// enum ImGuiHoveredFlags_
extern "C" const int const_ImGuiHoveredFlags_None;
extern "C" const int const_ImGuiHoveredFlags_ChildWindows;
extern "C" const int const_ImGuiHoveredFlags_RootWindow;
extern "C" const int const_ImGuiHoveredFlags_AnyWindow;
extern "C" const int const_ImGuiHoveredFlags_NoPopupHierarchy;
extern "C" const int const_ImGuiHoveredFlags_DockHierarchy;
extern "C" const int const_ImGuiHoveredFlags_AllowWhenBlockedByPopup;
extern "C" const int const_ImGuiHoveredFlags_AllowWhenBlockedByActiveItem;
extern "C" const int const_ImGuiHoveredFlags_AllowWhenOverlapped;
extern "C" const int const_ImGuiHoveredFlags_AllowWhenDisabled;
extern "C" const int const_ImGuiHoveredFlags_NoNavOverride;
extern "C" const int const_ImGuiHoveredFlags_RectOnly;
extern "C" const int const_ImGuiHoveredFlags_RootAndChildWindows;
// enum ImGuiDockNodeFlags_
extern "C" const int const_ImGuiDockNodeFlags_None;
extern "C" const int const_ImGuiDockNodeFlags_KeepAliveOnly;
extern "C" const int const_ImGuiDockNodeFlags_NoDockingInCentralNode;
extern "C" const int const_ImGuiDockNodeFlags_PassthruCentralNode;
extern "C" const int const_ImGuiDockNodeFlags_NoSplit;
extern "C" const int const_ImGuiDockNodeFlags_NoResize;
extern "C" const int const_ImGuiDockNodeFlags_AutoHideTabBar;
// enum ImGuiDragDropFlags_
extern "C" const int const_ImGuiDragDropFlags_None;
extern "C" const int const_ImGuiDragDropFlags_SourceNoPreviewTooltip;
extern "C" const int const_ImGuiDragDropFlags_SourceNoDisableHover;
extern "C" const int const_ImGuiDragDropFlags_SourceNoHoldToOpenOthers;
extern "C" const int const_ImGuiDragDropFlags_SourceAllowNullID;
extern "C" const int const_ImGuiDragDropFlags_SourceExtern;
extern "C" const int const_ImGuiDragDropFlags_SourceAutoExpirePayload;
extern "C" const int const_ImGuiDragDropFlags_AcceptBeforeDelivery;
extern "C" const int const_ImGuiDragDropFlags_AcceptNoDrawDefaultRect;
extern "C" const int const_ImGuiDragDropFlags_AcceptNoPreviewTooltip;
extern "C" const int const_ImGuiDragDropFlags_AcceptPeekOnly;
// enum ImGuiDataType_
extern "C" const int const_ImGuiDataType_S8;
extern "C" const int const_ImGuiDataType_U8;
extern "C" const int const_ImGuiDataType_S16;
extern "C" const int const_ImGuiDataType_U16;
extern "C" const int const_ImGuiDataType_S32;
extern "C" const int const_ImGuiDataType_U32;
extern "C" const int const_ImGuiDataType_S64;
extern "C" const int const_ImGuiDataType_U64;
extern "C" const int const_ImGuiDataType_Float;
extern "C" const int const_ImGuiDataType_Double;
extern "C" const int const_ImGuiDataType_COUNT;
// enum ImGuiDir_
extern "C" const int const_ImGuiDir_None;
extern "C" const int const_ImGuiDir_Left;
extern "C" const int const_ImGuiDir_Right;
extern "C" const int const_ImGuiDir_Up;
extern "C" const int const_ImGuiDir_Down;
extern "C" const int const_ImGuiDir_COUNT;
// enum ImGuiSortDirection_
extern "C" const int const_ImGuiSortDirection_None;
extern "C" const int const_ImGuiSortDirection_Ascending;
extern "C" const int const_ImGuiSortDirection_Descending;
// enum ImGuiKey_
extern "C" const int const_ImGuiKey_None;
extern "C" const int const_ImGuiKey_Tab;
extern "C" const int const_ImGuiKey_LeftArrow;
extern "C" const int const_ImGuiKey_RightArrow;
extern "C" const int const_ImGuiKey_UpArrow;
extern "C" const int const_ImGuiKey_DownArrow;
extern "C" const int const_ImGuiKey_PageUp;
extern "C" const int const_ImGuiKey_PageDown;
extern "C" const int const_ImGuiKey_Home;
extern "C" const int const_ImGuiKey_End;
extern "C" const int const_ImGuiKey_Insert;
extern "C" const int const_ImGuiKey_Delete;
extern "C" const int const_ImGuiKey_Backspace;
extern "C" const int const_ImGuiKey_Space;
extern "C" const int const_ImGuiKey_Enter;
extern "C" const int const_ImGuiKey_Escape;
extern "C" const int const_ImGuiKey_LeftCtrl;
extern "C" const int const_ImGuiKey_LeftShift;
extern "C" const int const_ImGuiKey_LeftAlt;
extern "C" const int const_ImGuiKey_LeftSuper;
extern "C" const int const_ImGuiKey_RightCtrl;
extern "C" const int const_ImGuiKey_RightShift;
extern "C" const int const_ImGuiKey_RightAlt;
extern "C" const int const_ImGuiKey_RightSuper;
extern "C" const int const_ImGuiKey_Menu;
extern "C" const int const_ImGuiKey_0;
extern "C" const int const_ImGuiKey_1;
extern "C" const int const_ImGuiKey_2;
extern "C" const int const_ImGuiKey_3;
extern "C" const int const_ImGuiKey_4;
extern "C" const int const_ImGuiKey_5;
extern "C" const int const_ImGuiKey_6;
extern "C" const int const_ImGuiKey_7;
extern "C" const int const_ImGuiKey_8;
extern "C" const int const_ImGuiKey_9;
extern "C" const int const_ImGuiKey_A;
extern "C" const int const_ImGuiKey_B;
extern "C" const int const_ImGuiKey_C;
extern "C" const int const_ImGuiKey_D;
extern "C" const int const_ImGuiKey_E;
extern "C" const int const_ImGuiKey_F;
extern "C" const int const_ImGuiKey_G;
extern "C" const int const_ImGuiKey_H;
extern "C" const int const_ImGuiKey_I;
extern "C" const int const_ImGuiKey_J;
extern "C" const int const_ImGuiKey_K;
extern "C" const int const_ImGuiKey_L;
extern "C" const int const_ImGuiKey_M;
extern "C" const int const_ImGuiKey_N;
extern "C" const int const_ImGuiKey_O;
extern "C" const int const_ImGuiKey_P;
extern "C" const int const_ImGuiKey_Q;
extern "C" const int const_ImGuiKey_R;
extern "C" const int const_ImGuiKey_S;
extern "C" const int const_ImGuiKey_T;
extern "C" const int const_ImGuiKey_U;
extern "C" const int const_ImGuiKey_V;
extern "C" const int const_ImGuiKey_W;
extern "C" const int const_ImGuiKey_X;
extern "C" const int const_ImGuiKey_Y;
extern "C" const int const_ImGuiKey_Z;
extern "C" const int const_ImGuiKey_F1;
extern "C" const int const_ImGuiKey_F2;
extern "C" const int const_ImGuiKey_F3;
extern "C" const int const_ImGuiKey_F4;
extern "C" const int const_ImGuiKey_F5;
extern "C" const int const_ImGuiKey_F6;
extern "C" const int const_ImGuiKey_F7;
extern "C" const int const_ImGuiKey_F8;
extern "C" const int const_ImGuiKey_F9;
extern "C" const int const_ImGuiKey_F10;
extern "C" const int const_ImGuiKey_F11;
extern "C" const int const_ImGuiKey_F12;
extern "C" const int const_ImGuiKey_Apostrophe;
extern "C" const int const_ImGuiKey_Comma;
extern "C" const int const_ImGuiKey_Minus;
extern "C" const int const_ImGuiKey_Period;
extern "C" const int const_ImGuiKey_Slash;
extern "C" const int const_ImGuiKey_Semicolon;
extern "C" const int const_ImGuiKey_Equal;
extern "C" const int const_ImGuiKey_LeftBracket;
extern "C" const int const_ImGuiKey_Backslash;
extern "C" const int const_ImGuiKey_RightBracket;
extern "C" const int const_ImGuiKey_GraveAccent;
extern "C" const int const_ImGuiKey_CapsLock;
extern "C" const int const_ImGuiKey_ScrollLock;
extern "C" const int const_ImGuiKey_NumLock;
extern "C" const int const_ImGuiKey_PrintScreen;
extern "C" const int const_ImGuiKey_Pause;
extern "C" const int const_ImGuiKey_Keypad0;
extern "C" const int const_ImGuiKey_Keypad1;
extern "C" const int const_ImGuiKey_Keypad2;
extern "C" const int const_ImGuiKey_Keypad3;
extern "C" const int const_ImGuiKey_Keypad4;
extern "C" const int const_ImGuiKey_Keypad5;
extern "C" const int const_ImGuiKey_Keypad6;
extern "C" const int const_ImGuiKey_Keypad7;
extern "C" const int const_ImGuiKey_Keypad8;
extern "C" const int const_ImGuiKey_Keypad9;
extern "C" const int const_ImGuiKey_KeypadDecimal;
extern "C" const int const_ImGuiKey_KeypadDivide;
extern "C" const int const_ImGuiKey_KeypadMultiply;
extern "C" const int const_ImGuiKey_KeypadSubtract;
extern "C" const int const_ImGuiKey_KeypadAdd;
extern "C" const int const_ImGuiKey_KeypadEnter;
extern "C" const int const_ImGuiKey_KeypadEqual;
extern "C" const int const_ImGuiKey_GamepadStart;
extern "C" const int const_ImGuiKey_GamepadBack;
extern "C" const int const_ImGuiKey_GamepadFaceUp;
extern "C" const int const_ImGuiKey_GamepadFaceDown;
extern "C" const int const_ImGuiKey_GamepadFaceLeft;
extern "C" const int const_ImGuiKey_GamepadFaceRight;
extern "C" const int const_ImGuiKey_GamepadDpadUp;
extern "C" const int const_ImGuiKey_GamepadDpadDown;
extern "C" const int const_ImGuiKey_GamepadDpadLeft;
extern "C" const int const_ImGuiKey_GamepadDpadRight;
extern "C" const int const_ImGuiKey_GamepadL1;
extern "C" const int const_ImGuiKey_GamepadR1;
extern "C" const int const_ImGuiKey_GamepadL2;
extern "C" const int const_ImGuiKey_GamepadR2;
extern "C" const int const_ImGuiKey_GamepadL3;
extern "C" const int const_ImGuiKey_GamepadR3;
extern "C" const int const_ImGuiKey_GamepadLStickUp;
extern "C" const int const_ImGuiKey_GamepadLStickDown;
extern "C" const int const_ImGuiKey_GamepadLStickLeft;
extern "C" const int const_ImGuiKey_GamepadLStickRight;
extern "C" const int const_ImGuiKey_GamepadRStickUp;
extern "C" const int const_ImGuiKey_GamepadRStickDown;
extern "C" const int const_ImGuiKey_GamepadRStickLeft;
extern "C" const int const_ImGuiKey_GamepadRStickRight;
extern "C" const int const_ImGuiKey_ModCtrl;
extern "C" const int const_ImGuiKey_ModShift;
extern "C" const int const_ImGuiKey_ModAlt;
extern "C" const int const_ImGuiKey_ModSuper;
// ImGuiKey_COUNT
extern "C" const int const_ImGuiKey_NamedKey_BEGIN;
extern "C" const int const_ImGuiKey_NamedKey_END;
// extern "C" const int const_ImGuiKey_NamedKey_COUNT;
// ImGuiKey_KeysData_SIZE
extern "C" const int const_ImGuiKey_KeysData_OFFSET;
// enum ImGuiModFlags_
extern "C" const int const_ImGuiModFlags_None;
extern "C" const int const_ImGuiModFlags_Ctrl;
extern "C" const int const_ImGuiModFlags_Shift;
extern "C" const int const_ImGuiModFlags_Alt;
extern "C" const int const_ImGuiModFlags_Super;
// enum ImGuiNavInput_
extern "C" const int const_ImGuiNavInput_Activate;
extern "C" const int const_ImGuiNavInput_Cancel;
extern "C" const int const_ImGuiNavInput_Input;
extern "C" const int const_ImGuiNavInput_Menu;
extern "C" const int const_ImGuiNavInput_DpadLeft;
extern "C" const int const_ImGuiNavInput_DpadRight;
extern "C" const int const_ImGuiNavInput_DpadUp;
extern "C" const int const_ImGuiNavInput_DpadDown;
extern "C" const int const_ImGuiNavInput_LStickLeft;
extern "C" const int const_ImGuiNavInput_LStickRight;
extern "C" const int const_ImGuiNavInput_LStickUp;
extern "C" const int const_ImGuiNavInput_LStickDown;
extern "C" const int const_ImGuiNavInput_FocusPrev;
extern "C" const int const_ImGuiNavInput_FocusNext;
extern "C" const int const_ImGuiNavInput_TweakSlow;
extern "C" const int const_ImGuiNavInput_TweakFast;
extern "C" const int const_ImGuiNavInput_KeyLeft_;
extern "C" const int const_ImGuiNavInput_KeyRight_;
extern "C" const int const_ImGuiNavInput_KeyUp_;
extern "C" const int const_ImGuiNavInput_KeyDown_;
// ImGuiNavInput_COUNT
// enum ImGuiConfigFlags_
extern "C" const int const_ImGuiConfigFlags_None;
extern "C" const int const_ImGuiConfigFlags_NavEnableKeyboard;
extern "C" const int const_ImGuiConfigFlags_NavEnableGamepad;
extern "C" const int const_ImGuiConfigFlags_NavEnableSetMousePos;
extern "C" const int const_ImGuiConfigFlags_NavNoCaptureKeyboard;
extern "C" const int const_ImGuiConfigFlags_NoMouse;
extern "C" const int const_ImGuiConfigFlags_NoMouseCursorChange;
extern "C" const int const_ImGuiConfigFlags_DockingEnable;
extern "C" const int const_ImGuiConfigFlags_ViewportsEnable;
extern "C" const int const_ImGuiConfigFlags_DpiEnableScaleViewports;
extern "C" const int const_ImGuiConfigFlags_DpiEnableScaleFonts;
extern "C" const int const_ImGuiConfigFlags_IsSRGB;
extern "C" const int const_ImGuiConfigFlags_IsTouchScreen;
// enum ImGuiBackendFlags_
extern "C" const int const_ImGuiBackendFlags_None;
extern "C" const int const_ImGuiBackendFlags_HasGamepad;
extern "C" const int const_ImGuiBackendFlags_HasMouseCursors;
extern "C" const int const_ImGuiBackendFlags_HasSetMousePos;
extern "C" const int const_ImGuiBackendFlags_RendererHasVtxOffset;
extern "C" const int const_ImGuiBackendFlags_PlatformHasViewports;
extern "C" const int const_ImGuiBackendFlags_HasMouseHoveredViewport;
extern "C" const int const_ImGuiBackendFlags_RendererHasViewports;
// enum ImGuiCol_
extern "C" const int const_ImGuiCol_Text;
extern "C" const int const_ImGuiCol_TextDisabled;
extern "C" const int const_ImGuiCol_WindowBg;
extern "C" const int const_ImGuiCol_ChildBg;
extern "C" const int const_ImGuiCol_PopupBg;
extern "C" const int const_ImGuiCol_Border;
extern "C" const int const_ImGuiCol_BorderShadow;
extern "C" const int const_ImGuiCol_FrameBg;
extern "C" const int const_ImGuiCol_FrameBgHovered;
extern "C" const int const_ImGuiCol_FrameBgActive;
extern "C" const int const_ImGuiCol_TitleBg;
extern "C" const int const_ImGuiCol_TitleBgActive;
extern "C" const int const_ImGuiCol_TitleBgCollapsed;
extern "C" const int const_ImGuiCol_MenuBarBg;
extern "C" const int const_ImGuiCol_ScrollbarBg;
extern "C" const int const_ImGuiCol_ScrollbarGrab;
extern "C" const int const_ImGuiCol_ScrollbarGrabHovered;
extern "C" const int const_ImGuiCol_ScrollbarGrabActive;
extern "C" const int const_ImGuiCol_CheckMark;
extern "C" const int const_ImGuiCol_SliderGrab;
extern "C" const int const_ImGuiCol_SliderGrabActive;
extern "C" const int const_ImGuiCol_Button;
extern "C" const int const_ImGuiCol_ButtonHovered;
extern "C" const int const_ImGuiCol_ButtonActive;
extern "C" const int const_ImGuiCol_Header;
extern "C" const int const_ImGuiCol_HeaderHovered;
extern "C" const int const_ImGuiCol_HeaderActive;
extern "C" const int const_ImGuiCol_Separator;
extern "C" const int const_ImGuiCol_SeparatorHovered;
extern "C" const int const_ImGuiCol_SeparatorActive;
extern "C" const int const_ImGuiCol_ResizeGrip;
extern "C" const int const_ImGuiCol_ResizeGripHovered;
extern "C" const int const_ImGuiCol_ResizeGripActive;
extern "C" const int const_ImGuiCol_Tab;
extern "C" const int const_ImGuiCol_TabHovered;
extern "C" const int const_ImGuiCol_TabActive;
extern "C" const int const_ImGuiCol_TabUnfocused;
extern "C" const int const_ImGuiCol_TabUnfocusedActive;
extern "C" const int const_ImGuiCol_DockingPreview;
extern "C" const int const_ImGuiCol_DockingEmptyBg;
extern "C" const int const_ImGuiCol_PlotLines;
extern "C" const int const_ImGuiCol_PlotLinesHovered;
extern "C" const int const_ImGuiCol_PlotHistogram;
extern "C" const int const_ImGuiCol_PlotHistogramHovered;
extern "C" const int const_ImGuiCol_TableHeaderBg;
extern "C" const int const_ImGuiCol_TableBorderStrong;
extern "C" const int const_ImGuiCol_TableBorderLight;
extern "C" const int const_ImGuiCol_TableRowBg;
extern "C" const int const_ImGuiCol_TableRowBgAlt;
extern "C" const int const_ImGuiCol_TextSelectedBg;
extern "C" const int const_ImGuiCol_DragDropTarget;
extern "C" const int const_ImGuiCol_NavHighlight;
extern "C" const int const_ImGuiCol_NavWindowingHighlight;
extern "C" const int const_ImGuiCol_NavWindowingDimBg;
extern "C" const int const_ImGuiCol_ModalWindowDimBg;
// extern "C" const int const_ImGuiCol_COUNT;
// enum ImGuiStyleVar_
extern "C" const int const_ImGuiStyleVar_Alpha;
extern "C" const int const_ImGuiStyleVar_DisabledAlpha;
extern "C" const int const_ImGuiStyleVar_WindowPadding;
extern "C" const int const_ImGuiStyleVar_WindowRounding;
extern "C" const int const_ImGuiStyleVar_WindowBorderSize;
extern "C" const int const_ImGuiStyleVar_WindowMinSize;
extern "C" const int const_ImGuiStyleVar_WindowTitleAlign;
extern "C" const int const_ImGuiStyleVar_ChildRounding;
extern "C" const int const_ImGuiStyleVar_ChildBorderSize;
extern "C" const int const_ImGuiStyleVar_PopupRounding;
extern "C" const int const_ImGuiStyleVar_PopupBorderSize;
extern "C" const int const_ImGuiStyleVar_FramePadding;
extern "C" const int const_ImGuiStyleVar_FrameRounding;
extern "C" const int const_ImGuiStyleVar_FrameBorderSize;
extern "C" const int const_ImGuiStyleVar_ItemSpacing;
extern "C" const int const_ImGuiStyleVar_ItemInnerSpacing;
extern "C" const int const_ImGuiStyleVar_IndentSpacing;
extern "C" const int const_ImGuiStyleVar_CellPadding;
extern "C" const int const_ImGuiStyleVar_ScrollbarSize;
extern "C" const int const_ImGuiStyleVar_ScrollbarRounding;
extern "C" const int const_ImGuiStyleVar_GrabMinSize;
extern "C" const int const_ImGuiStyleVar_GrabRounding;
extern "C" const int const_ImGuiStyleVar_TabRounding;
extern "C" const int const_ImGuiStyleVar_ButtonTextAlign;
extern "C" const int const_ImGuiStyleVar_SelectableTextAlign;
extern "C" const int const_ImGuiStyleVar_COUNT;
// enum ImGuiButtonFlags_
extern "C" const int const_ImGuiButtonFlags_None;
extern "C" const int const_ImGuiButtonFlags_MouseButtonLeft;
extern "C" const int const_ImGuiButtonFlags_MouseButtonRight;
extern "C" const int const_ImGuiButtonFlags_MouseButtonMiddle;
extern "C" const int const_ImGuiButtonFlags_MouseButtonMask_;
extern "C" const int const_ImGuiButtonFlags_MouseButtonDefault_;
// enum ImGuiColorEditFlags_
extern "C" const int const_ImGuiColorEditFlags_None;
extern "C" const int const_ImGuiColorEditFlags_NoAlpha;
extern "C" const int const_ImGuiColorEditFlags_NoPicker;
extern "C" const int const_ImGuiColorEditFlags_NoOptions;
extern "C" const int const_ImGuiColorEditFlags_NoSmallPreview;
extern "C" const int const_ImGuiColorEditFlags_NoInputs;
extern "C" const int const_ImGuiColorEditFlags_NoTooltip;
extern "C" const int const_ImGuiColorEditFlags_NoLabel;
extern "C" const int const_ImGuiColorEditFlags_NoSidePreview;
extern "C" const int const_ImGuiColorEditFlags_NoDragDrop;
extern "C" const int const_ImGuiColorEditFlags_NoBorder;
extern "C" const int const_ImGuiColorEditFlags_AlphaBar;
extern "C" const int const_ImGuiColorEditFlags_AlphaPreview;
extern "C" const int const_ImGuiColorEditFlags_AlphaPreviewHalf;
extern "C" const int const_ImGuiColorEditFlags_HDR;
extern "C" const int const_ImGuiColorEditFlags_DisplayRGB;
extern "C" const int const_ImGuiColorEditFlags_DisplayHSV;
extern "C" const int const_ImGuiColorEditFlags_DisplayHex;
extern "C" const int const_ImGuiColorEditFlags_Uint8;
extern "C" const int const_ImGuiColorEditFlags_Float;
extern "C" const int const_ImGuiColorEditFlags_PickerHueBar;
extern "C" const int const_ImGuiColorEditFlags_PickerHueWheel;
extern "C" const int const_ImGuiColorEditFlags_InputRGB;
extern "C" const int const_ImGuiColorEditFlags_InputHSV;
extern "C" const int const_ImGuiColorEditFlags_DefaultOptions_;
extern "C" const int const_ImGuiColorEditFlags_DisplayMask_;
extern "C" const int const_ImGuiColorEditFlags_DataTypeMask_;
extern "C" const int const_ImGuiColorEditFlags_PickerMask_;
extern "C" const int const_ImGuiColorEditFlags_InputMask_;
// enum ImGuiSliderFlags_
extern "C" const int const_ImGuiSliderFlags_None;
extern "C" const int const_ImGuiSliderFlags_AlwaysClamp;
extern "C" const int const_ImGuiSliderFlags_Logarithmic;
extern "C" const int const_ImGuiSliderFlags_NoRoundToFormat;
extern "C" const int const_ImGuiSliderFlags_NoInput;
extern "C" const int const_ImGuiSliderFlags_InvalidMask_;
// enum ImGuiMouseButton_
extern "C" const int const_ImGuiMouseButton_Left;
extern "C" const int const_ImGuiMouseButton_Right;
extern "C" const int const_ImGuiMouseButton_Middle;
extern "C" const int const_ImGuiMouseButton_COUNT;
// enum ImGuiMouseCursor_
extern "C" const int const_ImGuiMouseCursor_None;
extern "C" const int const_ImGuiMouseCursor_Arrow;
extern "C" const int const_ImGuiMouseCursor_TextInput;
extern "C" const int const_ImGuiMouseCursor_ResizeAll;
extern "C" const int const_ImGuiMouseCursor_ResizeNS;
extern "C" const int const_ImGuiMouseCursor_ResizeEW;
extern "C" const int const_ImGuiMouseCursor_ResizeNESW;
extern "C" const int const_ImGuiMouseCursor_ResizeNWSE;
extern "C" const int const_ImGuiMouseCursor_Hand;
extern "C" const int const_ImGuiMouseCursor_NotAllowed;
extern "C" const int const_ImGuiMouseCursor_COUNT;
// enum ImGuiCond_
extern "C" const int const_ImGuiCond_None;
extern "C" const int const_ImGuiCond_Always;
extern "C" const int const_ImGuiCond_Once;
extern "C" const int const_ImGuiCond_FirstUseEver;
extern "C" const int const_ImGuiCond_Appearing;
// enum ImDrawFlags_
extern "C" const int const_ImDrawFlags_None;
extern "C" const int const_ImDrawFlags_Closed;
extern "C" const int const_ImDrawFlags_RoundCornersTopLeft;
extern "C" const int const_ImDrawFlags_RoundCornersTopRight;
extern "C" const int const_ImDrawFlags_RoundCornersBottomLeft;
extern "C" const int const_ImDrawFlags_RoundCornersBottomRight;
extern "C" const int const_ImDrawFlags_RoundCornersNone;
extern "C" const int const_ImDrawFlags_RoundCornersTop;
extern "C" const int const_ImDrawFlags_RoundCornersBottom;
extern "C" const int const_ImDrawFlags_RoundCornersLeft;
extern "C" const int const_ImDrawFlags_RoundCornersRight;
extern "C" const int const_ImDrawFlags_RoundCornersAll;
extern "C" const int const_ImDrawFlags_RoundCornersDefault_;
extern "C" const int const_ImDrawFlags_RoundCornersMask_;
// enum ImDrawListFlags_
extern "C" const int const_ImDrawListFlags_None;
extern "C" const int const_ImDrawListFlags_AntiAliasedLines;
extern "C" const int const_ImDrawListFlags_AntiAliasedLinesUseTex;
extern "C" const int const_ImDrawListFlags_AntiAliasedFill;
extern "C" const int const_ImDrawListFlags_AllowVtxOffset;
// enum ImFontAtlasFlags_
extern "C" const int const_ImFontAtlasFlags_None;
extern "C" const int const_ImFontAtlasFlags_NoPowerOfTwoHeight;
extern "C" const int const_ImFontAtlasFlags_NoMouseCursors;
extern "C" const int const_ImFontAtlasFlags_NoBakedLines;
// enum ImGuiViewportFlags_
extern "C" const int const_ImGuiViewportFlags_None;
extern "C" const int const_ImGuiViewportFlags_IsPlatformWindow;
extern "C" const int const_ImGuiViewportFlags_IsPlatformMonitor;
extern "C" const int const_ImGuiViewportFlags_OwnedByApp;
extern "C" const int const_ImGuiViewportFlags_NoDecoration;
extern "C" const int const_ImGuiViewportFlags_NoTaskBarIcon;
extern "C" const int const_ImGuiViewportFlags_NoFocusOnAppearing;
extern "C" const int const_ImGuiViewportFlags_NoFocusOnClick;
extern "C" const int const_ImGuiViewportFlags_NoInputs;
extern "C" const int const_ImGuiViewportFlags_NoRendererClear;
extern "C" const int const_ImGuiViewportFlags_TopMost;
extern "C" const int const_ImGuiViewportFlags_Minimized;
extern "C" const int const_ImGuiViewportFlags_NoAutoMerge;
extern "C" const int const_ImGuiViewportFlags_CanHostOtherWindows;
// enum ImGuiItemFlags_
extern "C" const int const_ImGuiItemFlags_None;
extern "C" const int const_ImGuiItemFlags_NoTabStop;
extern "C" const int const_ImGuiItemFlags_ButtonRepeat;
extern "C" const int const_ImGuiItemFlags_Disabled;
extern "C" const int const_ImGuiItemFlags_NoNav;
extern "C" const int const_ImGuiItemFlags_NoNavDefaultFocus;
extern "C" const int const_ImGuiItemFlags_SelectableDontClosePopup;
extern "C" const int const_ImGuiItemFlags_MixedValue;
extern "C" const int const_ImGuiItemFlags_ReadOnly;
extern "C" const int const_ImGuiItemFlags_Inputable;
// enum ImGuiItemStatusFlags_
extern "C" const int const_ImGuiItemStatusFlags_None;
extern "C" const int const_ImGuiItemStatusFlags_HoveredRect;
extern "C" const int const_ImGuiItemStatusFlags_HasDisplayRect;
extern "C" const int const_ImGuiItemStatusFlags_Edited;
extern "C" const int const_ImGuiItemStatusFlags_ToggledSelection;
extern "C" const int const_ImGuiItemStatusFlags_ToggledOpen;
extern "C" const int const_ImGuiItemStatusFlags_HasDeactivated;
extern "C" const int const_ImGuiItemStatusFlags_Deactivated;
extern "C" const int const_ImGuiItemStatusFlags_HoveredWindow;
extern "C" const int const_ImGuiItemStatusFlags_FocusedByTabbing;
// enum ImGuiInputTextFlagsPrivate_
extern "C" const int const_ImGuiInputTextFlags_Multiline;
extern "C" const int const_ImGuiInputTextFlags_NoMarkEdited;
extern "C" const int const_ImGuiInputTextFlags_MergedItem;
// enum ImGuiButtonFlagsPrivate_
extern "C" const int const_ImGuiButtonFlags_PressedOnClick;
extern "C" const int const_ImGuiButtonFlags_PressedOnClickRelease;
extern "C" const int const_ImGuiButtonFlags_PressedOnClickReleaseAnywhere;
extern "C" const int const_ImGuiButtonFlags_PressedOnRelease;
extern "C" const int const_ImGuiButtonFlags_PressedOnDoubleClick;
extern "C" const int const_ImGuiButtonFlags_PressedOnDragDropHold;
extern "C" const int const_ImGuiButtonFlags_Repeat;
extern "C" const int const_ImGuiButtonFlags_FlattenChildren;
extern "C" const int const_ImGuiButtonFlags_AllowItemOverlap;
extern "C" const int const_ImGuiButtonFlags_DontClosePopups;
extern "C" const int const_ImGuiButtonFlags_AlignTextBaseLine;
extern "C" const int const_ImGuiButtonFlags_NoKeyModifiers;
extern "C" const int const_ImGuiButtonFlags_NoHoldingActiveId;
extern "C" const int const_ImGuiButtonFlags_NoNavFocus;
extern "C" const int const_ImGuiButtonFlags_NoHoveredOnFocus;
extern "C" const int const_ImGuiButtonFlags_PressedOnMask_;
extern "C" const int const_ImGuiButtonFlags_PressedOnDefault_;
// enum ImGuiComboFlagsPrivate_
extern "C" const int const_ImGuiComboFlags_CustomPreview;
// enum ImGuiSliderFlagsPrivate_
extern "C" const int const_ImGuiSliderFlags_Vertical;
extern "C" const int const_ImGuiSliderFlags_ReadOnly;
// enum ImGuiSelectableFlagsPrivate_
extern "C" const int const_ImGuiSelectableFlags_NoHoldingActiveID;
extern "C" const int const_ImGuiSelectableFlags_SelectOnNav;
extern "C" const int const_ImGuiSelectableFlags_SelectOnClick;
extern "C" const int const_ImGuiSelectableFlags_SelectOnRelease;
extern "C" const int const_ImGuiSelectableFlags_SpanAvailWidth;
extern "C" const int const_ImGuiSelectableFlags_DrawHoveredWhenHeld;
extern "C" const int const_ImGuiSelectableFlags_SetNavIdOnHover;
extern "C" const int const_ImGuiSelectableFlags_NoPadWithHalfSpacing;
// enum ImGuiTreeNodeFlagsPrivate_
extern "C" const int const_ImGuiTreeNodeFlags_ClipLabelForTrailingButton;
// enum ImGuiSeparatorFlags_
extern "C" const int const_ImGuiSeparatorFlags_None;
extern "C" const int const_ImGuiSeparatorFlags_Horizontal;
extern "C" const int const_ImGuiSeparatorFlags_Vertical;
extern "C" const int const_ImGuiSeparatorFlags_SpanAllColumns;
// enum ImGuiTextFlags_
extern "C" const int const_ImGuiTextFlags_None;
extern "C" const int const_ImGuiTextFlags_NoWidthForLargeClippedText;
// enum ImGuiTooltipFlags_
extern "C" const int const_ImGuiTooltipFlags_None;
extern "C" const int const_ImGuiTooltipFlags_OverridePreviousTooltip;
// enum ImGuiLayoutType_
extern "C" const int const_ImGuiLayoutType_Horizontal;
extern "C" const int const_ImGuiLayoutType_Vertical;
// enum ImGuiLogType
extern "C" const int const_ImGuiLogType_None;
extern "C" const int const_ImGuiLogType_TTY;
extern "C" const int const_ImGuiLogType_File;
extern "C" const int const_ImGuiLogType_Buffer;
extern "C" const int const_ImGuiLogType_Clipboard;
// enum ImGuiAxis
extern "C" const int const_ImGuiAxis_None;
extern "C" const int const_ImGuiAxis_X;
extern "C" const int const_ImGuiAxis_Y;
// enum ImGuiPlotType
extern "C" const int const_ImGuiPlotType_Lines;
extern "C" const int const_ImGuiPlotType_Histogram;
// enum ImGuiPopupPositionPolicy
extern "C" const int const_ImGuiPopupPositionPolicy_Default;
extern "C" const int const_ImGuiPopupPositionPolicy_ComboBox;
extern "C" const int const_ImGuiPopupPositionPolicy_Tooltip;
// enum ImGuiDataTypePrivate_
extern "C" const int const_ImGuiDataType_String;
extern "C" const int const_ImGuiDataType_Pointer;
extern "C" const int const_ImGuiDataType_ID;
// enum ImGuiNextWindowDataFlags_
extern "C" const int const_ImGuiNextWindowDataFlags_None;
extern "C" const int const_ImGuiNextWindowDataFlags_HasPos;
extern "C" const int const_ImGuiNextWindowDataFlags_HasSize;
extern "C" const int const_ImGuiNextWindowDataFlags_HasContentSize;
extern "C" const int const_ImGuiNextWindowDataFlags_HasCollapsed;
extern "C" const int const_ImGuiNextWindowDataFlags_HasSizeConstraint;
extern "C" const int const_ImGuiNextWindowDataFlags_HasFocus;
extern "C" const int const_ImGuiNextWindowDataFlags_HasBgAlpha;
extern "C" const int const_ImGuiNextWindowDataFlags_HasScroll;
extern "C" const int const_ImGuiNextWindowDataFlags_HasViewport;
extern "C" const int const_ImGuiNextWindowDataFlags_HasDock;
extern "C" const int const_ImGuiNextWindowDataFlags_HasWindowClass;
// enum ImGuiNextItemDataFlags_
extern "C" const int const_ImGuiNextItemDataFlags_None;
extern "C" const int const_ImGuiNextItemDataFlags_HasWidth;
extern "C" const int const_ImGuiNextItemDataFlags_HasOpen;
// enum ImGuiKeyPrivate_
extern "C" const int const_ImGuiKey_LegacyNativeKey_BEGIN;
extern "C" const int const_ImGuiKey_LegacyNativeKey_END;
extern "C" const int const_ImGuiKey_Gamepad_BEGIN;
extern "C" const int const_ImGuiKey_Gamepad_END;
// enum ImGuiInputEventType
extern "C" const int const_ImGuiInputEventType_None;
extern "C" const int const_ImGuiInputEventType_MousePos;
extern "C" const int const_ImGuiInputEventType_MouseWheel;
extern "C" const int const_ImGuiInputEventType_MouseButton;
extern "C" const int const_ImGuiInputEventType_MouseViewport;
extern "C" const int const_ImGuiInputEventType_Key;
extern "C" const int const_ImGuiInputEventType_Text;
extern "C" const int const_ImGuiInputEventType_Focus;
extern "C" const int const_ImGuiInputEventType_COUNT;
// enum ImGuiInputSource
extern "C" const int const_ImGuiInputSource_None;
extern "C" const int const_ImGuiInputSource_Mouse;
extern "C" const int const_ImGuiInputSource_Keyboard;
extern "C" const int const_ImGuiInputSource_Gamepad;
extern "C" const int const_ImGuiInputSource_Clipboard;
extern "C" const int const_ImGuiInputSource_Nav;
extern "C" const int const_ImGuiInputSource_COUNT;
// enum ImGuiNavReadMode
extern "C" const int const_ImGuiNavReadMode_Down;
extern "C" const int const_ImGuiNavReadMode_Pressed;
extern "C" const int const_ImGuiNavReadMode_Released;
extern "C" const int const_ImGuiNavReadMode_Repeat;
extern "C" const int const_ImGuiNavReadMode_RepeatSlow;
extern "C" const int const_ImGuiNavReadMode_RepeatFast;
// enum ImGuiActivateFlags_
extern "C" const int const_ImGuiActivateFlags_None;
extern "C" const int const_ImGuiActivateFlags_PreferInput;
extern "C" const int const_ImGuiActivateFlags_PreferTweak;
extern "C" const int const_ImGuiActivateFlags_TryToPreserveState;
// enum ImGuiScrollFlags_
extern "C" const int const_ImGuiScrollFlags_None;
extern "C" const int const_ImGuiScrollFlags_KeepVisibleEdgeX;
extern "C" const int const_ImGuiScrollFlags_KeepVisibleEdgeY;
extern "C" const int const_ImGuiScrollFlags_KeepVisibleCenterX;
extern "C" const int const_ImGuiScrollFlags_KeepVisibleCenterY;
extern "C" const int const_ImGuiScrollFlags_AlwaysCenterX;
extern "C" const int const_ImGuiScrollFlags_AlwaysCenterY;
extern "C" const int const_ImGuiScrollFlags_NoScrollParent;
extern "C" const int const_ImGuiScrollFlags_MaskX_;
extern "C" const int const_ImGuiScrollFlags_MaskY_;
// enum ImGuiNavHighlightFlags_
extern "C" const int const_ImGuiNavHighlightFlags_None;
extern "C" const int const_ImGuiNavHighlightFlags_TypeDefault;
extern "C" const int const_ImGuiNavHighlightFlags_TypeThin;
extern "C" const int const_ImGuiNavHighlightFlags_AlwaysDraw;
extern "C" const int const_ImGuiNavHighlightFlags_NoRounding;
// enum ImGuiNavDirSourceFlags_
extern "C" const int const_ImGuiNavDirSourceFlags_None;
extern "C" const int const_ImGuiNavDirSourceFlags_RawKeyboard;
extern "C" const int const_ImGuiNavDirSourceFlags_Keyboard;
extern "C" const int const_ImGuiNavDirSourceFlags_PadDPad;
extern "C" const int const_ImGuiNavDirSourceFlags_PadLStick;
// enum ImGuiNavMoveFlags_
extern "C" const int const_ImGuiNavMoveFlags_None;
extern "C" const int const_ImGuiNavMoveFlags_LoopX;
extern "C" const int const_ImGuiNavMoveFlags_LoopY;
extern "C" const int const_ImGuiNavMoveFlags_WrapX;
extern "C" const int const_ImGuiNavMoveFlags_WrapY;
extern "C" const int const_ImGuiNavMoveFlags_AllowCurrentNavId;
extern "C" const int const_ImGuiNavMoveFlags_AlsoScoreVisibleSet;
extern "C" const int const_ImGuiNavMoveFlags_ScrollToEdgeY;
extern "C" const int const_ImGuiNavMoveFlags_Forwarded;
extern "C" const int const_ImGuiNavMoveFlags_DebugNoResult;
extern "C" const int const_ImGuiNavMoveFlags_FocusApi;
extern "C" const int const_ImGuiNavMoveFlags_Tabbing;
extern "C" const int const_ImGuiNavMoveFlags_Activate;
extern "C" const int const_ImGuiNavMoveFlags_DontSetNavHighlight;
// enum ImGuiNavLayer
extern "C" const int const_ImGuiNavLayer_Main;
extern "C" const int const_ImGuiNavLayer_Menu;
// extern "C" const int const_ImGuiNavLayer_COUNT;
// enum ImGuiOldColumnFlags_
extern "C" const int const_ImGuiOldColumnFlags_None;
extern "C" const int const_ImGuiOldColumnFlags_NoBorder;
extern "C" const int const_ImGuiOldColumnFlags_NoResize;
extern "C" const int const_ImGuiOldColumnFlags_NoPreserveWidths;
extern "C" const int const_ImGuiOldColumnFlags_NoForceWithinWindow;
extern "C" const int const_ImGuiOldColumnFlags_GrowParentContentsSize;
// enum ImGuiDockNodeFlagsPrivate_
extern "C" const int const_ImGuiDockNodeFlags_DockSpace;
extern "C" const int const_ImGuiDockNodeFlags_CentralNode;
extern "C" const int const_ImGuiDockNodeFlags_NoTabBar;
extern "C" const int const_ImGuiDockNodeFlags_HiddenTabBar;
extern "C" const int const_ImGuiDockNodeFlags_NoWindowMenuButton;
extern "C" const int const_ImGuiDockNodeFlags_NoCloseButton;
extern "C" const int const_ImGuiDockNodeFlags_NoDocking;
extern "C" const int const_ImGuiDockNodeFlags_NoDockingSplitMe;
extern "C" const int const_ImGuiDockNodeFlags_NoDockingSplitOther;
extern "C" const int const_ImGuiDockNodeFlags_NoDockingOverMe;
extern "C" const int const_ImGuiDockNodeFlags_NoDockingOverOther;
extern "C" const int const_ImGuiDockNodeFlags_NoDockingOverEmpty;
extern "C" const int const_ImGuiDockNodeFlags_NoResizeX;
extern "C" const int const_ImGuiDockNodeFlags_NoResizeY;
extern "C" const int const_ImGuiDockNodeFlags_SharedFlagsInheritMask_;
extern "C" const int const_ImGuiDockNodeFlags_NoResizeFlagsMask_;
extern "C" const int const_ImGuiDockNodeFlags_LocalFlagsMask_;
extern "C" const int const_ImGuiDockNodeFlags_LocalFlagsTransferMask_;
extern "C" const int const_ImGuiDockNodeFlags_SavedFlagsMask_;
// enum ImGuiDataAuthority_
extern "C" const int const_ImGuiDataAuthority_Auto;
extern "C" const int const_ImGuiDataAuthority_DockNode;
extern "C" const int const_ImGuiDataAuthority_Window;
// enum ImGuiDockNodeState
extern "C" const int const_ImGuiDockNodeState_Unknown;
extern "C" const int const_ImGuiDockNodeState_HostWindowHiddenBecauseSingleWindow;
extern "C" const int const_ImGuiDockNodeState_HostWindowHiddenBecauseWindowsAreResizing;
extern "C" const int const_ImGuiDockNodeState_HostWindowVisible;
// enum ImGuiWindowDockStyleCol
extern "C" const int const_ImGuiWindowDockStyleCol_Text;
extern "C" const int const_ImGuiWindowDockStyleCol_Tab;
extern "C" const int const_ImGuiWindowDockStyleCol_TabHovered;
extern "C" const int const_ImGuiWindowDockStyleCol_TabActive;
extern "C" const int const_ImGuiWindowDockStyleCol_TabUnfocused;
extern "C" const int const_ImGuiWindowDockStyleCol_TabUnfocusedActive;
// extern "C" const int const_ImGuiWindowDockStyleCol_COUNT;
// enum ImGuiDebugLogFlags_
extern "C" const int const_ImGuiDebugLogFlags_None;
extern "C" const int const_ImGuiDebugLogFlags_EventActiveId;
extern "C" const int const_ImGuiDebugLogFlags_EventFocus;
extern "C" const int const_ImGuiDebugLogFlags_EventPopup;
extern "C" const int const_ImGuiDebugLogFlags_EventNav;
extern "C" const int const_ImGuiDebugLogFlags_EventIO;
extern "C" const int const_ImGuiDebugLogFlags_EventDocking;
extern "C" const int const_ImGuiDebugLogFlags_EventViewport;
extern "C" const int const_ImGuiDebugLogFlags_EventMask_;
extern "C" const int const_ImGuiDebugLogFlags_OutputToTTY;
// enum ImGuiContextHookType
extern "C" const int const_ImGuiContextHookType_NewFramePre;
extern "C" const int const_ImGuiContextHookType_NewFramePost;
extern "C" const int const_ImGuiContextHookType_EndFramePre;
extern "C" const int const_ImGuiContextHookType_EndFramePost;
extern "C" const int const_ImGuiContextHookType_RenderPre;
extern "C" const int const_ImGuiContextHookType_RenderPost;
extern "C" const int const_ImGuiContextHookType_Shutdown;
extern "C" const int const_ImGuiContextHookType_PendingRemoval_;
// enum ImGuiTabBarFlagsPrivate_
extern "C" const int const_ImGuiTabBarFlags_DockNode;
extern "C" const int const_ImGuiTabBarFlags_IsFocused;
extern "C" const int const_ImGuiTabBarFlags_SaveSettings;
// enum ImGuiTabItemFlagsPrivate_
extern "C" const int const_ImGuiTabItemFlags_SectionMask_;
extern "C" const int const_ImGuiTabItemFlags_NoCloseButton;
extern "C" const int const_ImGuiTabItemFlags_Button;
extern "C" const int const_ImGuiTabItemFlags_Unsorted;
extern "C" const int const_ImGuiTabItemFlags_Preview;
// enum ImGuiFreeTypeBuilderFlags
extern "C" const int const_ImGuiFreeTypeBuilderFlags_NoHinting;
extern "C" const int const_ImGuiFreeTypeBuilderFlags_NoAutoHint;
extern "C" const int const_ImGuiFreeTypeBuilderFlags_ForceAutoHint;
extern "C" const int const_ImGuiFreeTypeBuilderFlags_LightHinting;
extern "C" const int const_ImGuiFreeTypeBuilderFlags_MonoHinting;
extern "C" const int const_ImGuiFreeTypeBuilderFlags_Bold;
extern "C" const int const_ImGuiFreeTypeBuilderFlags_Oblique;
extern "C" const int const_ImGuiFreeTypeBuilderFlags_Monochrome;
extern "C" const int const_ImGuiFreeTypeBuilderFlags_LoadColor;
extern "C" const int const_ImGuiFreeTypeBuilderFlags_Bitmap;

//// filedialog
// enum ImGuiFileDialogFlags;
extern "C" const int const_ImGuiFileDialogFlags_None;
extern "C" const int const_ImGuiFileDialogFlags_ConfirmOverwrite;
extern "C" const int const_ImGuiFileDialogFlags_DontShowHiddenFiles;
extern "C" const int const_ImGuiFileDialogFlags_DisableCreateDirectoryButton;
extern "C" const int const_ImGuiFileDialogFlags_HideColumnType;
extern "C" const int const_ImGuiFileDialogFlags_HideColumnSize;
extern "C" const int const_ImGuiFileDialogFlags_HideColumnDate;
extern "C" const int const_ImGuiFileDialogFlags_NoDialog;
extern "C" const int const_ImGuiFileDialogFlags_ReadOnlyFileNameField;
extern "C" const int const_ImGuiFileDialogFlags_CaseInsensitiveExtention;
extern "C" const int const_ImGuiFileDialogFlags_Modal;
// enum IGFD_FileStyle
extern "C" const int const_IGFD_FileStyle_None;
extern "C" const int const_IGFD_FileStyleByTypeFile;
extern "C" const int const_IGFD_FileStyleByTypeDir;
extern "C" const int const_IGFD_FileStyleByTypeLink;
extern "C" const int const_IGFD_FileStyleByExtention;
extern "C" const int const_IGFD_FileStyleByFullName;
extern "C" const int const_IGFD_FileStyleByContainedInFullName;

//// implot
// enum ImAxis_
extern "C" const int const_ImAxis_X1;
extern "C" const int const_ImAxis_X2;
extern "C" const int const_ImAxis_X3;
extern "C" const int const_ImAxis_Y1;
extern "C" const int const_ImAxis_Y2;
extern "C" const int const_ImAxis_Y3;
// enum ImPlotFlags_
extern "C" const int const_ImPlotFlags_None;
extern "C" const int const_ImPlotFlags_NoTitle;
extern "C" const int const_ImPlotFlags_NoLegend;
extern "C" const int const_ImPlotFlags_NoMouseText;
extern "C" const int const_ImPlotFlags_NoInputs;
extern "C" const int const_ImPlotFlags_NoMenus;
extern "C" const int const_ImPlotFlags_NoBoxSelect;
extern "C" const int const_ImPlotFlags_NoChild;
extern "C" const int const_ImPlotFlags_NoFrame;
extern "C" const int const_ImPlotFlags_Equal;
extern "C" const int const_ImPlotFlags_Crosshairs;
extern "C" const int const_ImPlotFlags_CanvasOnly;
// enum ImPlotAxisFlags_
extern "C" const int const_ImPlotAxisFlags_None;
extern "C" const int const_ImPlotAxisFlags_NoLabel;
extern "C" const int const_ImPlotAxisFlags_NoGridLines;
extern "C" const int const_ImPlotAxisFlags_NoTickMarks;
extern "C" const int const_ImPlotAxisFlags_NoTickLabels;
extern "C" const int const_ImPlotAxisFlags_NoInitialFit;
extern "C" const int const_ImPlotAxisFlags_NoMenus;
extern "C" const int const_ImPlotAxisFlags_NoSideSwitch;
extern "C" const int const_ImPlotAxisFlags_NoHighlight;
extern "C" const int const_ImPlotAxisFlags_Opposite;
extern "C" const int const_ImPlotAxisFlags_Foreground;
extern "C" const int const_ImPlotAxisFlags_Invert;
extern "C" const int const_ImPlotAxisFlags_AutoFit;
extern "C" const int const_ImPlotAxisFlags_RangeFit;
extern "C" const int const_ImPlotAxisFlags_PanStretch;
extern "C" const int const_ImPlotAxisFlags_LockMin;
extern "C" const int const_ImPlotAxisFlags_LockMax;
extern "C" const int const_ImPlotAxisFlags_Lock;
extern "C" const int const_ImPlotAxisFlags_NoDecorations;
extern "C" const int const_ImPlotAxisFlags_AuxDefault;
// enum ImPlotLineFlags_
extern "C" const int const_ImPlotLineFlags_None;
extern "C" const int const_ImPlotLineFlags_Segments;
extern "C" const int const_ImPlotLineFlags_Loop;
extern "C" const int const_ImPlotLineFlags_SkipNaN;
extern "C" const int const_ImPlotLineFlags_NoClip;
extern "C" const int const_ImPlotLineFlags_Shaded;
// enum ImPlotMarker_
extern "C" const int const_ImPlotMarker_None;
extern "C" const int const_ImPlotMarker_Circle;
extern "C" const int const_ImPlotMarker_Square;
extern "C" const int const_ImPlotMarker_Diamond;
extern "C" const int const_ImPlotMarker_Up;
extern "C" const int const_ImPlotMarker_Down;
extern "C" const int const_ImPlotMarker_Left;
extern "C" const int const_ImPlotMarker_Right;
extern "C" const int const_ImPlotMarker_Cross;
extern "C" const int const_ImPlotMarker_Plus;
extern "C" const int const_ImPlotMarker_Asterisk;
// const int const_ImPlotMarker_COUNT

//// The fonts for the critic2 GUI
extern "C" const char *const_font_dejavu_base85_ptr;

// The definitions in glcorearb.h
#define GL_DEPTH_BUFFER_BIT               0x00000100
#define GL_STENCIL_BUFFER_BIT             0x00000400
#define GL_COLOR_BUFFER_BIT               0x00004000
#define GL_FALSE                          0
#define GL_TRUE                           1
#define GL_POINTS                         0x0000
#define GL_LINES                          0x0001
#define GL_LINE_LOOP                      0x0002
#define GL_LINE_STRIP                     0x0003
#define GL_TRIANGLES                      0x0004
#define GL_TRIANGLE_STRIP                 0x0005
#define GL_TRIANGLE_FAN                   0x0006
#define GL_QUADS                          0x0007
#define GL_NEVER                          0x0200
#define GL_LESS                           0x0201
#define GL_EQUAL                          0x0202
#define GL_LEQUAL                         0x0203
#define GL_GREATER                        0x0204
#define GL_NOTEQUAL                       0x0205
#define GL_GEQUAL                         0x0206
#define GL_ALWAYS                         0x0207
#define GL_ZERO                           0
#define GL_ONE                            1
#define GL_SRC_COLOR                      0x0300
#define GL_ONE_MINUS_SRC_COLOR            0x0301
#define GL_SRC_ALPHA                      0x0302
#define GL_ONE_MINUS_SRC_ALPHA            0x0303
#define GL_DST_ALPHA                      0x0304
#define GL_ONE_MINUS_DST_ALPHA            0x0305
#define GL_DST_COLOR                      0x0306
#define GL_ONE_MINUS_DST_COLOR            0x0307
#define GL_SRC_ALPHA_SATURATE             0x0308
#define GL_NONE                           0
#define GL_FRONT_LEFT                     0x0400
#define GL_FRONT_RIGHT                    0x0401
#define GL_BACK_LEFT                      0x0402
#define GL_BACK_RIGHT                     0x0403
#define GL_FRONT                          0x0404
#define GL_BACK                           0x0405
#define GL_LEFT                           0x0406
#define GL_RIGHT                          0x0407
#define GL_FRONT_AND_BACK                 0x0408
#define GL_NO_ERROR                       0
#define GL_INVALID_ENUM                   0x0500
#define GL_INVALID_VALUE                  0x0501
#define GL_INVALID_OPERATION              0x0502
#define GL_OUT_OF_MEMORY                  0x0505
#define GL_CW                             0x0900
#define GL_CCW                            0x0901
#define GL_POINT_SIZE                     0x0B11
#define GL_POINT_SIZE_RANGE               0x0B12
#define GL_POINT_SIZE_GRANULARITY         0x0B13
#define GL_LINE_SMOOTH                    0x0B20
#define GL_LINE_WIDTH                     0x0B21
#define GL_LINE_WIDTH_RANGE               0x0B22
#define GL_LINE_WIDTH_GRANULARITY         0x0B23
#define GL_POLYGON_MODE                   0x0B40
#define GL_POLYGON_SMOOTH                 0x0B41
#define GL_CULL_FACE                      0x0B44
#define GL_CULL_FACE_MODE                 0x0B45
#define GL_FRONT_FACE                     0x0B46
#define GL_DEPTH_RANGE                    0x0B70
#define GL_DEPTH_TEST                     0x0B71
#define GL_DEPTH_WRITEMASK                0x0B72
#define GL_DEPTH_CLEAR_VALUE              0x0B73
#define GL_DEPTH_FUNC                     0x0B74
#define GL_STENCIL_TEST                   0x0B90
#define GL_STENCIL_CLEAR_VALUE            0x0B91
#define GL_STENCIL_FUNC                   0x0B92
#define GL_STENCIL_VALUE_MASK             0x0B93
#define GL_STENCIL_FAIL                   0x0B94
#define GL_STENCIL_PASS_DEPTH_FAIL        0x0B95
#define GL_STENCIL_PASS_DEPTH_PASS        0x0B96
#define GL_STENCIL_REF                    0x0B97
#define GL_STENCIL_WRITEMASK              0x0B98
#define GL_VIEWPORT                       0x0BA2
#define GL_DITHER                         0x0BD0
#define GL_BLEND_DST                      0x0BE0
#define GL_BLEND_SRC                      0x0BE1
#define GL_BLEND                          0x0BE2
#define GL_LOGIC_OP_MODE                  0x0BF0
#define GL_COLOR_LOGIC_OP                 0x0BF2
#define GL_DRAW_BUFFER                    0x0C01
#define GL_READ_BUFFER                    0x0C02
#define GL_SCISSOR_BOX                    0x0C10
#define GL_SCISSOR_TEST                   0x0C11
#define GL_COLOR_CLEAR_VALUE              0x0C22
#define GL_COLOR_WRITEMASK                0x0C23
#define GL_DOUBLEBUFFER                   0x0C32
#define GL_STEREO                         0x0C33
#define GL_LINE_SMOOTH_HINT               0x0C52
#define GL_POLYGON_SMOOTH_HINT            0x0C53
#define GL_UNPACK_SWAP_BYTES              0x0CF0
#define GL_UNPACK_LSB_FIRST               0x0CF1
#define GL_UNPACK_ROW_LENGTH              0x0CF2
#define GL_UNPACK_SKIP_ROWS               0x0CF3
#define GL_UNPACK_SKIP_PIXELS             0x0CF4
#define GL_UNPACK_ALIGNMENT               0x0CF5
#define GL_PACK_SWAP_BYTES                0x0D00
#define GL_PACK_LSB_FIRST                 0x0D01
#define GL_PACK_ROW_LENGTH                0x0D02
#define GL_PACK_SKIP_ROWS                 0x0D03
#define GL_PACK_SKIP_PIXELS               0x0D04
#define GL_PACK_ALIGNMENT                 0x0D05
#define GL_MAX_TEXTURE_SIZE               0x0D33
#define GL_MAX_VIEWPORT_DIMS              0x0D3A
#define GL_SUBPIXEL_BITS                  0x0D50
#define GL_TEXTURE_1D                     0x0DE0
#define GL_TEXTURE_2D                     0x0DE1
#define GL_POLYGON_OFFSET_UNITS           0x2A00
#define GL_POLYGON_OFFSET_POINT           0x2A01
#define GL_POLYGON_OFFSET_LINE            0x2A02
#define GL_POLYGON_OFFSET_FILL            0x8037
#define GL_POLYGON_OFFSET_FACTOR          0x8038
#define GL_TEXTURE_BINDING_1D             0x8068
#define GL_TEXTURE_BINDING_2D             0x8069
#define GL_TEXTURE_WIDTH                  0x1000
#define GL_TEXTURE_HEIGHT                 0x1001
#define GL_TEXTURE_INTERNAL_FORMAT        0x1003
#define GL_TEXTURE_BORDER_COLOR           0x1004
#define GL_TEXTURE_RED_SIZE               0x805C
#define GL_TEXTURE_GREEN_SIZE             0x805D
#define GL_TEXTURE_BLUE_SIZE              0x805E
#define GL_TEXTURE_ALPHA_SIZE             0x805F
#define GL_DONT_CARE                      0x1100
#define GL_FASTEST                        0x1101
#define GL_NICEST                         0x1102
#define GL_BYTE                           0x1400
#define GL_UNSIGNED_BYTE                  0x1401
#define GL_SHORT                          0x1402
#define GL_UNSIGNED_SHORT                 0x1403
#define GL_INT                            0x1404
#define GL_UNSIGNED_INT                   0x1405
#define GL_FLOAT                          0x1406
#define GL_DOUBLE                         0x140A
#define GL_STACK_OVERFLOW                 0x0503
#define GL_STACK_UNDERFLOW                0x0504
#define GL_CLEAR                          0x1500
#define GL_AND                            0x1501
#define GL_AND_REVERSE                    0x1502
#define GL_COPY                           0x1503
#define GL_AND_INVERTED                   0x1504
#define GL_NOOP                           0x1505
#define GL_XOR                            0x1506
#define GL_OR                             0x1507
#define GL_NOR                            0x1508
#define GL_EQUIV                          0x1509
#define GL_INVERT                         0x150A
#define GL_OR_REVERSE                     0x150B
#define GL_COPY_INVERTED                  0x150C
#define GL_OR_INVERTED                    0x150D
#define GL_NAND                           0x150E
#define GL_SET                            0x150F
#define GL_TEXTURE                        0x1702
#define GL_COLOR                          0x1800
#define GL_DEPTH                          0x1801
#define GL_STENCIL                        0x1802
#define GL_STENCIL_INDEX                  0x1901
#define GL_DEPTH_COMPONENT                0x1902
#define GL_RED                            0x1903
#define GL_GREEN                          0x1904
#define GL_BLUE                           0x1905
#define GL_ALPHA                          0x1906
#define GL_RGB                            0x1907
#define GL_RGBA                           0x1908
#define GL_POINT                          0x1B00
#define GL_LINE                           0x1B01
#define GL_FILL                           0x1B02
#define GL_KEEP                           0x1E00
#define GL_REPLACE                        0x1E01
#define GL_INCR                           0x1E02
#define GL_DECR                           0x1E03
#define GL_VENDOR                         0x1F00
#define GL_RENDERER                       0x1F01
#define GL_VERSION                        0x1F02
#define GL_EXTENSIONS                     0x1F03
#define GL_NEAREST                        0x2600
#define GL_LINEAR                         0x2601
#define GL_NEAREST_MIPMAP_NEAREST         0x2700
#define GL_LINEAR_MIPMAP_NEAREST          0x2701
#define GL_NEAREST_MIPMAP_LINEAR          0x2702
#define GL_LINEAR_MIPMAP_LINEAR           0x2703
#define GL_TEXTURE_MAG_FILTER             0x2800
#define GL_TEXTURE_MIN_FILTER             0x2801
#define GL_TEXTURE_WRAP_S                 0x2802
#define GL_TEXTURE_WRAP_T                 0x2803
#define GL_PROXY_TEXTURE_1D               0x8063
#define GL_PROXY_TEXTURE_2D               0x8064
#define GL_REPEAT                         0x2901
#define GL_R3_G3_B2                       0x2A10
#define GL_RGB4                           0x804F
#define GL_RGB5                           0x8050
#define GL_RGB8                           0x8051
#define GL_RGB10                          0x8052
#define GL_RGB12                          0x8053
#define GL_RGB16                          0x8054
#define GL_RGBA2                          0x8055
#define GL_RGBA4                          0x8056
#define GL_RGB5_A1                        0x8057
#define GL_RGBA8                          0x8058
#define GL_RGB10_A2                       0x8059
#define GL_RGBA12                         0x805A
#define GL_RGBA16                         0x805B
#define GL_UNSIGNED_BYTE_3_3_2            0x8032
#define GL_UNSIGNED_SHORT_4_4_4_4         0x8033
#define GL_UNSIGNED_SHORT_5_5_5_1         0x8034
#define GL_UNSIGNED_INT_8_8_8_8           0x8035
#define GL_UNSIGNED_INT_10_10_10_2        0x8036
#define GL_TEXTURE_BINDING_3D             0x806A
#define GL_PACK_SKIP_IMAGES               0x806B
#define GL_PACK_IMAGE_HEIGHT              0x806C
#define GL_UNPACK_SKIP_IMAGES             0x806D
#define GL_UNPACK_IMAGE_HEIGHT            0x806E
#define GL_TEXTURE_3D                     0x806F
#define GL_PROXY_TEXTURE_3D               0x8070
#define GL_TEXTURE_DEPTH                  0x8071
#define GL_TEXTURE_WRAP_R                 0x8072
#define GL_MAX_3D_TEXTURE_SIZE            0x8073
#define GL_UNSIGNED_BYTE_2_3_3_REV        0x8362
#define GL_UNSIGNED_SHORT_5_6_5           0x8363
#define GL_UNSIGNED_SHORT_5_6_5_REV       0x8364
#define GL_UNSIGNED_SHORT_4_4_4_4_REV     0x8365
#define GL_UNSIGNED_SHORT_1_5_5_5_REV     0x8366
#define GL_UNSIGNED_INT_8_8_8_8_REV       0x8367
#define GL_UNSIGNED_INT_2_10_10_10_REV    0x8368
#define GL_BGR                            0x80E0
#define GL_BGRA                           0x80E1
#define GL_MAX_ELEMENTS_VERTICES          0x80E8
#define GL_MAX_ELEMENTS_INDICES           0x80E9
#define GL_CLAMP_TO_EDGE                  0x812F
#define GL_TEXTURE_MIN_LOD                0x813A
#define GL_TEXTURE_MAX_LOD                0x813B
#define GL_TEXTURE_BASE_LEVEL             0x813C
#define GL_TEXTURE_MAX_LEVEL              0x813D
#define GL_SMOOTH_POINT_SIZE_RANGE        0x0B12
#define GL_SMOOTH_POINT_SIZE_GRANULARITY  0x0B13
#define GL_SMOOTH_LINE_WIDTH_RANGE        0x0B22
#define GL_SMOOTH_LINE_WIDTH_GRANULARITY  0x0B23
#define GL_ALIASED_LINE_WIDTH_RANGE       0x846E
#define GL_CONSTANT_COLOR                 0x8001
#define GL_ONE_MINUS_CONSTANT_COLOR       0x8002
#define GL_CONSTANT_ALPHA                 0x8003
#define GL_ONE_MINUS_CONSTANT_ALPHA       0x8004
#define GL_BLEND_COLOR                    0x8005
#define GL_FUNC_ADD                       0x8006
#define GL_MIN                            0x8007
#define GL_MAX                            0x8008
#define GL_BLEND_EQUATION                 0x8009
#define GL_FUNC_SUBTRACT                  0x800A
#define GL_FUNC_REVERSE_SUBTRACT          0x800B
#define GL_TEXTURE0                       0x84C0
#define GL_TEXTURE1                       0x84C1
#define GL_TEXTURE2                       0x84C2
#define GL_TEXTURE3                       0x84C3
#define GL_TEXTURE4                       0x84C4
#define GL_TEXTURE5                       0x84C5
#define GL_TEXTURE6                       0x84C6
#define GL_TEXTURE7                       0x84C7
#define GL_TEXTURE8                       0x84C8
#define GL_TEXTURE9                       0x84C9
#define GL_TEXTURE10                      0x84CA
#define GL_TEXTURE11                      0x84CB
#define GL_TEXTURE12                      0x84CC
#define GL_TEXTURE13                      0x84CD
#define GL_TEXTURE14                      0x84CE
#define GL_TEXTURE15                      0x84CF
#define GL_TEXTURE16                      0x84D0
#define GL_TEXTURE17                      0x84D1
#define GL_TEXTURE18                      0x84D2
#define GL_TEXTURE19                      0x84D3
#define GL_TEXTURE20                      0x84D4
#define GL_TEXTURE21                      0x84D5
#define GL_TEXTURE22                      0x84D6
#define GL_TEXTURE23                      0x84D7
#define GL_TEXTURE24                      0x84D8
#define GL_TEXTURE25                      0x84D9
#define GL_TEXTURE26                      0x84DA
#define GL_TEXTURE27                      0x84DB
#define GL_TEXTURE28                      0x84DC
#define GL_TEXTURE29                      0x84DD
#define GL_TEXTURE30                      0x84DE
#define GL_TEXTURE31                      0x84DF
#define GL_ACTIVE_TEXTURE                 0x84E0
#define GL_MULTISAMPLE                    0x809D
#define GL_SAMPLE_ALPHA_TO_COVERAGE       0x809E
#define GL_SAMPLE_ALPHA_TO_ONE            0x809F
#define GL_SAMPLE_COVERAGE                0x80A0
#define GL_SAMPLE_BUFFERS                 0x80A8
#define GL_SAMPLES                        0x80A9
#define GL_SAMPLE_COVERAGE_VALUE          0x80AA
#define GL_SAMPLE_COVERAGE_INVERT         0x80AB
#define GL_TEXTURE_CUBE_MAP               0x8513
#define GL_TEXTURE_BINDING_CUBE_MAP       0x8514
#define GL_TEXTURE_CUBE_MAP_POSITIVE_X    0x8515
#define GL_TEXTURE_CUBE_MAP_NEGATIVE_X    0x8516
#define GL_TEXTURE_CUBE_MAP_POSITIVE_Y    0x8517
#define GL_TEXTURE_CUBE_MAP_NEGATIVE_Y    0x8518
#define GL_TEXTURE_CUBE_MAP_POSITIVE_Z    0x8519
#define GL_TEXTURE_CUBE_MAP_NEGATIVE_Z    0x851A
#define GL_PROXY_TEXTURE_CUBE_MAP         0x851B
#define GL_MAX_CUBE_MAP_TEXTURE_SIZE      0x851C
#define GL_COMPRESSED_RGB                 0x84ED
#define GL_COMPRESSED_RGBA                0x84EE
#define GL_TEXTURE_COMPRESSION_HINT       0x84EF
#define GL_TEXTURE_COMPRESSED_IMAGE_SIZE  0x86A0
#define GL_TEXTURE_COMPRESSED             0x86A1
#define GL_NUM_COMPRESSED_TEXTURE_FORMATS 0x86A2
#define GL_COMPRESSED_TEXTURE_FORMATS     0x86A3
#define GL_CLAMP_TO_BORDER                0x812D
#define GL_BLEND_DST_RGB                  0x80C8
#define GL_BLEND_SRC_RGB                  0x80C9
#define GL_BLEND_DST_ALPHA                0x80CA
#define GL_BLEND_SRC_ALPHA                0x80CB
#define GL_POINT_FADE_THRESHOLD_SIZE      0x8128
#define GL_DEPTH_COMPONENT16              0x81A5
#define GL_DEPTH_COMPONENT24              0x81A6
#define GL_DEPTH_COMPONENT32              0x81A7
#define GL_MIRRORED_REPEAT                0x8370
#define GL_MAX_TEXTURE_LOD_BIAS           0x84FD
#define GL_TEXTURE_LOD_BIAS               0x8501
#define GL_INCR_WRAP                      0x8507
#define GL_DECR_WRAP                      0x8508
#define GL_TEXTURE_DEPTH_SIZE             0x884A
#define GL_TEXTURE_COMPARE_MODE           0x884C
#define GL_TEXTURE_COMPARE_FUNC           0x884D
#define GL_BUFFER_SIZE                    0x8764
#define GL_BUFFER_USAGE                   0x8765
#define GL_QUERY_COUNTER_BITS             0x8864
#define GL_CURRENT_QUERY                  0x8865
#define GL_QUERY_RESULT                   0x8866
#define GL_QUERY_RESULT_AVAILABLE         0x8867
#define GL_ARRAY_BUFFER                   0x8892
#define GL_ELEMENT_ARRAY_BUFFER           0x8893
#define GL_ARRAY_BUFFER_BINDING           0x8894
#define GL_ELEMENT_ARRAY_BUFFER_BINDING   0x8895
#define GL_VERTEX_ATTRIB_ARRAY_BUFFER_BINDING 0x889F
#define GL_READ_ONLY                      0x88B8
#define GL_WRITE_ONLY                     0x88B9
#define GL_READ_WRITE                     0x88BA
#define GL_BUFFER_ACCESS                  0x88BB
#define GL_BUFFER_MAPPED                  0x88BC
#define GL_BUFFER_MAP_POINTER             0x88BD
#define GL_STREAM_DRAW                    0x88E0
#define GL_STREAM_READ                    0x88E1
#define GL_STREAM_COPY                    0x88E2
#define GL_STATIC_DRAW                    0x88E4
#define GL_STATIC_READ                    0x88E5
#define GL_STATIC_COPY                    0x88E6
#define GL_DYNAMIC_DRAW                   0x88E8
#define GL_DYNAMIC_READ                   0x88E9
#define GL_DYNAMIC_COPY                   0x88EA
#define GL_SAMPLES_PASSED                 0x8914
#define GL_SRC1_ALPHA                     0x8589
#define GL_BLEND_EQUATION_RGB             0x8009
#define GL_VERTEX_ATTRIB_ARRAY_ENABLED    0x8622
#define GL_VERTEX_ATTRIB_ARRAY_SIZE       0x8623
#define GL_VERTEX_ATTRIB_ARRAY_STRIDE     0x8624
#define GL_VERTEX_ATTRIB_ARRAY_TYPE       0x8625
#define GL_CURRENT_VERTEX_ATTRIB          0x8626
#define GL_VERTEX_PROGRAM_POINT_SIZE      0x8642
#define GL_VERTEX_ATTRIB_ARRAY_POINTER    0x8645
#define GL_STENCIL_BACK_FUNC              0x8800
#define GL_STENCIL_BACK_FAIL              0x8801
#define GL_STENCIL_BACK_PASS_DEPTH_FAIL   0x8802
#define GL_STENCIL_BACK_PASS_DEPTH_PASS   0x8803
#define GL_MAX_DRAW_BUFFERS               0x8824
#define GL_DRAW_BUFFER0                   0x8825
#define GL_DRAW_BUFFER1                   0x8826
#define GL_DRAW_BUFFER2                   0x8827
#define GL_DRAW_BUFFER3                   0x8828
#define GL_DRAW_BUFFER4                   0x8829
#define GL_DRAW_BUFFER5                   0x882A
#define GL_DRAW_BUFFER6                   0x882B
#define GL_DRAW_BUFFER7                   0x882C
#define GL_DRAW_BUFFER8                   0x882D
#define GL_DRAW_BUFFER9                   0x882E
#define GL_DRAW_BUFFER10                  0x882F
#define GL_DRAW_BUFFER11                  0x8830
#define GL_DRAW_BUFFER12                  0x8831
#define GL_DRAW_BUFFER13                  0x8832
#define GL_DRAW_BUFFER14                  0x8833
#define GL_DRAW_BUFFER15                  0x8834
#define GL_BLEND_EQUATION_ALPHA           0x883D
#define GL_MAX_VERTEX_ATTRIBS             0x8869
#define GL_VERTEX_ATTRIB_ARRAY_NORMALIZED 0x886A
#define GL_MAX_TEXTURE_IMAGE_UNITS        0x8872
#define GL_FRAGMENT_SHADER                0x8B30
#define GL_VERTEX_SHADER                  0x8B31
#define GL_MAX_FRAGMENT_UNIFORM_COMPONENTS 0x8B49
#define GL_MAX_VERTEX_UNIFORM_COMPONENTS  0x8B4A
#define GL_MAX_VARYING_FLOATS             0x8B4B
#define GL_MAX_VERTEX_TEXTURE_IMAGE_UNITS 0x8B4C
#define GL_MAX_COMBINED_TEXTURE_IMAGE_UNITS 0x8B4D
#define GL_SHADER_TYPE                    0x8B4F
#define GL_FLOAT_VEC2                     0x8B50
#define GL_FLOAT_VEC3                     0x8B51
#define GL_FLOAT_VEC4                     0x8B52
#define GL_INT_VEC2                       0x8B53
#define GL_INT_VEC3                       0x8B54
#define GL_INT_VEC4                       0x8B55
#define GL_BOOL                           0x8B56
#define GL_BOOL_VEC2                      0x8B57
#define GL_BOOL_VEC3                      0x8B58
#define GL_BOOL_VEC4                      0x8B59
#define GL_FLOAT_MAT2                     0x8B5A
#define GL_FLOAT_MAT3                     0x8B5B
#define GL_FLOAT_MAT4                     0x8B5C
#define GL_SAMPLER_1D                     0x8B5D
#define GL_SAMPLER_2D                     0x8B5E
#define GL_SAMPLER_3D                     0x8B5F
#define GL_SAMPLER_CUBE                   0x8B60
#define GL_SAMPLER_1D_SHADOW              0x8B61
#define GL_SAMPLER_2D_SHADOW              0x8B62
#define GL_DELETE_STATUS                  0x8B80
#define GL_COMPILE_STATUS                 0x8B81
#define GL_LINK_STATUS                    0x8B82
#define GL_VALIDATE_STATUS                0x8B83
#define GL_INFO_LOG_LENGTH                0x8B84
#define GL_ATTACHED_SHADERS               0x8B85
#define GL_ACTIVE_UNIFORMS                0x8B86
#define GL_ACTIVE_UNIFORM_MAX_LENGTH      0x8B87
#define GL_SHADER_SOURCE_LENGTH           0x8B88
#define GL_ACTIVE_ATTRIBUTES              0x8B89
#define GL_ACTIVE_ATTRIBUTE_MAX_LENGTH    0x8B8A
#define GL_FRAGMENT_SHADER_DERIVATIVE_HINT 0x8B8B
#define GL_SHADING_LANGUAGE_VERSION       0x8B8C
#define GL_CURRENT_PROGRAM                0x8B8D
#define GL_POINT_SPRITE_COORD_ORIGIN      0x8CA0
#define GL_LOWER_LEFT                     0x8CA1
#define GL_UPPER_LEFT                     0x8CA2
#define GL_STENCIL_BACK_REF               0x8CA3
#define GL_STENCIL_BACK_VALUE_MASK        0x8CA4
#define GL_STENCIL_BACK_WRITEMASK         0x8CA5
#define GL_PIXEL_PACK_BUFFER              0x88EB
#define GL_PIXEL_UNPACK_BUFFER            0x88EC
#define GL_PIXEL_PACK_BUFFER_BINDING      0x88ED
#define GL_PIXEL_UNPACK_BUFFER_BINDING    0x88EF
#define GL_FLOAT_MAT2x3                   0x8B65
#define GL_FLOAT_MAT2x4                   0x8B66
#define GL_FLOAT_MAT3x2                   0x8B67
#define GL_FLOAT_MAT3x4                   0x8B68
#define GL_FLOAT_MAT4x2                   0x8B69
#define GL_FLOAT_MAT4x3                   0x8B6A
#define GL_SRGB                           0x8C40
#define GL_SRGB8                          0x8C41
#define GL_SRGB_ALPHA                     0x8C42
#define GL_SRGB8_ALPHA8                   0x8C43
#define GL_COMPRESSED_SRGB                0x8C48
#define GL_COMPRESSED_SRGB_ALPHA          0x8C49
#define GL_COMPARE_REF_TO_TEXTURE         0x884E
#define GL_CLIP_DISTANCE0                 0x3000
#define GL_CLIP_DISTANCE1                 0x3001
#define GL_CLIP_DISTANCE2                 0x3002
#define GL_CLIP_DISTANCE3                 0x3003
#define GL_CLIP_DISTANCE4                 0x3004
#define GL_CLIP_DISTANCE5                 0x3005
#define GL_CLIP_DISTANCE6                 0x3006
#define GL_CLIP_DISTANCE7                 0x3007
#define GL_MAX_CLIP_DISTANCES             0x0D32
#define GL_MAJOR_VERSION                  0x821B
#define GL_MINOR_VERSION                  0x821C
#define GL_NUM_EXTENSIONS                 0x821D
#define GL_CONTEXT_FLAGS                  0x821E
#define GL_COMPRESSED_RED                 0x8225
#define GL_COMPRESSED_RG                  0x8226
#define GL_CONTEXT_FLAG_FORWARD_COMPATIBLE_BIT 0x0001
#define GL_RGBA32F                        0x8814
#define GL_RGB32F                         0x8815
#define GL_RGBA16F                        0x881A
#define GL_RGB16F                         0x881B
#define GL_VERTEX_ATTRIB_ARRAY_INTEGER    0x88FD
#define GL_MAX_ARRAY_TEXTURE_LAYERS       0x88FF
#define GL_MIN_PROGRAM_TEXEL_OFFSET       0x8904
#define GL_MAX_PROGRAM_TEXEL_OFFSET       0x8905
#define GL_CLAMP_READ_COLOR               0x891C
#define GL_FIXED_ONLY                     0x891D
#define GL_MAX_VARYING_COMPONENTS         0x8B4B
#define GL_TEXTURE_1D_ARRAY               0x8C18
#define GL_PROXY_TEXTURE_1D_ARRAY         0x8C19
#define GL_TEXTURE_2D_ARRAY               0x8C1A
#define GL_PROXY_TEXTURE_2D_ARRAY         0x8C1B
#define GL_TEXTURE_BINDING_1D_ARRAY       0x8C1C
#define GL_TEXTURE_BINDING_2D_ARRAY       0x8C1D
#define GL_R11F_G11F_B10F                 0x8C3A
#define GL_UNSIGNED_INT_10F_11F_11F_REV   0x8C3B
#define GL_RGB9_E5                        0x8C3D
#define GL_UNSIGNED_INT_5_9_9_9_REV       0x8C3E
#define GL_TEXTURE_SHARED_SIZE            0x8C3F
#define GL_TRANSFORM_FEEDBACK_VARYING_MAX_LENGTH 0x8C76
#define GL_TRANSFORM_FEEDBACK_BUFFER_MODE 0x8C7F
#define GL_MAX_TRANSFORM_FEEDBACK_SEPARATE_COMPONENTS 0x8C80
#define GL_TRANSFORM_FEEDBACK_VARYINGS    0x8C83
#define GL_TRANSFORM_FEEDBACK_BUFFER_START 0x8C84
#define GL_TRANSFORM_FEEDBACK_BUFFER_SIZE 0x8C85
#define GL_PRIMITIVES_GENERATED           0x8C87
#define GL_TRANSFORM_FEEDBACK_PRIMITIVES_WRITTEN 0x8C88
#define GL_RASTERIZER_DISCARD             0x8C89
#define GL_MAX_TRANSFORM_FEEDBACK_INTERLEAVED_COMPONENTS 0x8C8A
#define GL_MAX_TRANSFORM_FEEDBACK_SEPARATE_ATTRIBS 0x8C8B
#define GL_INTERLEAVED_ATTRIBS            0x8C8C
#define GL_SEPARATE_ATTRIBS               0x8C8D
#define GL_TRANSFORM_FEEDBACK_BUFFER      0x8C8E
#define GL_TRANSFORM_FEEDBACK_BUFFER_BINDING 0x8C8F
#define GL_RGBA32UI                       0x8D70
#define GL_RGB32UI                        0x8D71
#define GL_RGBA16UI                       0x8D76
#define GL_RGB16UI                        0x8D77
#define GL_RGBA8UI                        0x8D7C
#define GL_RGB8UI                         0x8D7D
#define GL_RGBA32I                        0x8D82
#define GL_RGB32I                         0x8D83
#define GL_RGBA16I                        0x8D88
#define GL_RGB16I                         0x8D89
#define GL_RGBA8I                         0x8D8E
#define GL_RGB8I                          0x8D8F
#define GL_RED_INTEGER                    0x8D94
#define GL_GREEN_INTEGER                  0x8D95
#define GL_BLUE_INTEGER                   0x8D96
#define GL_RGB_INTEGER                    0x8D98
#define GL_RGBA_INTEGER                   0x8D99
#define GL_BGR_INTEGER                    0x8D9A
#define GL_BGRA_INTEGER                   0x8D9B
#define GL_SAMPLER_1D_ARRAY               0x8DC0
#define GL_SAMPLER_2D_ARRAY               0x8DC1
#define GL_SAMPLER_1D_ARRAY_SHADOW        0x8DC3
#define GL_SAMPLER_2D_ARRAY_SHADOW        0x8DC4
#define GL_SAMPLER_CUBE_SHADOW            0x8DC5
#define GL_UNSIGNED_INT_VEC2              0x8DC6
#define GL_UNSIGNED_INT_VEC3              0x8DC7
#define GL_UNSIGNED_INT_VEC4              0x8DC8
#define GL_INT_SAMPLER_1D                 0x8DC9
#define GL_INT_SAMPLER_2D                 0x8DCA
#define GL_INT_SAMPLER_3D                 0x8DCB
#define GL_INT_SAMPLER_CUBE               0x8DCC
#define GL_INT_SAMPLER_1D_ARRAY           0x8DCE
#define GL_INT_SAMPLER_2D_ARRAY           0x8DCF
#define GL_UNSIGNED_INT_SAMPLER_1D        0x8DD1
#define GL_UNSIGNED_INT_SAMPLER_2D        0x8DD2
#define GL_UNSIGNED_INT_SAMPLER_3D        0x8DD3
#define GL_UNSIGNED_INT_SAMPLER_CUBE      0x8DD4
#define GL_UNSIGNED_INT_SAMPLER_1D_ARRAY  0x8DD6
#define GL_UNSIGNED_INT_SAMPLER_2D_ARRAY  0x8DD7
#define GL_QUERY_WAIT                     0x8E13
#define GL_QUERY_NO_WAIT                  0x8E14
#define GL_QUERY_BY_REGION_WAIT           0x8E15
#define GL_QUERY_BY_REGION_NO_WAIT        0x8E16
#define GL_BUFFER_ACCESS_FLAGS            0x911F
#define GL_BUFFER_MAP_LENGTH              0x9120
#define GL_BUFFER_MAP_OFFSET              0x9121
#define GL_SAMPLER_2D_RECT                0x8B63
#define GL_SAMPLER_2D_RECT_SHADOW         0x8B64
#define GL_SAMPLER_BUFFER                 0x8DC2
#define GL_INT_SAMPLER_2D_RECT            0x8DCD
#define GL_INT_SAMPLER_BUFFER             0x8DD0
#define GL_UNSIGNED_INT_SAMPLER_2D_RECT   0x8DD5
#define GL_UNSIGNED_INT_SAMPLER_BUFFER    0x8DD8
#define GL_TEXTURE_BUFFER                 0x8C2A
#define GL_MAX_TEXTURE_BUFFER_SIZE        0x8C2B
#define GL_TEXTURE_BINDING_BUFFER         0x8C2C
#define GL_TEXTURE_BUFFER_DATA_STORE_BINDING 0x8C2D
#define GL_TEXTURE_RECTANGLE              0x84F5
#define GL_TEXTURE_BINDING_RECTANGLE      0x84F6
#define GL_PROXY_TEXTURE_RECTANGLE        0x84F7
#define GL_MAX_RECTANGLE_TEXTURE_SIZE     0x84F8
#define GL_RED_SNORM                      0x8F90
#define GL_RG_SNORM                       0x8F91
#define GL_RGB_SNORM                      0x8F92
#define GL_RGBA_SNORM                     0x8F93
#define GL_R8_SNORM                       0x8F94
#define GL_RG8_SNORM                      0x8F95
#define GL_RGB8_SNORM                     0x8F96
#define GL_RGBA8_SNORM                    0x8F97
#define GL_R16_SNORM                      0x8F98
#define GL_RG16_SNORM                     0x8F99
#define GL_RGB16_SNORM                    0x8F9A
#define GL_RGBA16_SNORM                   0x8F9B
#define GL_SIGNED_NORMALIZED              0x8F9C
#define GL_PRIMITIVE_RESTART              0x8F9D
#define GL_PRIMITIVE_RESTART_INDEX        0x8F9E
#define GL_CONTEXT_CORE_PROFILE_BIT       0x00000001
#define GL_CONTEXT_COMPATIBILITY_PROFILE_BIT 0x00000002
#define GL_LINES_ADJACENCY                0x000A
#define GL_LINE_STRIP_ADJACENCY           0x000B
#define GL_TRIANGLES_ADJACENCY            0x000C
#define GL_TRIANGLE_STRIP_ADJACENCY       0x000D
#define GL_PROGRAM_POINT_SIZE             0x8642
#define GL_MAX_GEOMETRY_TEXTURE_IMAGE_UNITS 0x8C29
#define GL_FRAMEBUFFER_ATTACHMENT_LAYERED 0x8DA7
#define GL_FRAMEBUFFER_INCOMPLETE_LAYER_TARGETS 0x8DA8
#define GL_GEOMETRY_SHADER                0x8DD9
#define GL_GEOMETRY_VERTICES_OUT          0x8916
#define GL_GEOMETRY_INPUT_TYPE            0x8917
#define GL_GEOMETRY_OUTPUT_TYPE           0x8918
#define GL_MAX_GEOMETRY_UNIFORM_COMPONENTS 0x8DDF
#define GL_MAX_GEOMETRY_OUTPUT_VERTICES   0x8DE0
#define GL_MAX_GEOMETRY_TOTAL_OUTPUT_COMPONENTS 0x8DE1
#define GL_MAX_VERTEX_OUTPUT_COMPONENTS   0x9122
#define GL_MAX_GEOMETRY_INPUT_COMPONENTS  0x9123
#define GL_MAX_GEOMETRY_OUTPUT_COMPONENTS 0x9124
#define GL_MAX_FRAGMENT_INPUT_COMPONENTS  0x9125
#define GL_CONTEXT_PROFILE_MASK           0x9126
#define GL_VERTEX_ATTRIB_ARRAY_DIVISOR    0x88FE
#define GL_SAMPLE_SHADING                 0x8C36
#define GL_MIN_SAMPLE_SHADING_VALUE       0x8C37
#define GL_MIN_PROGRAM_TEXTURE_GATHER_OFFSET 0x8E5E
#define GL_MAX_PROGRAM_TEXTURE_GATHER_OFFSET 0x8E5F
#define GL_TEXTURE_CUBE_MAP_ARRAY         0x9009
#define GL_TEXTURE_BINDING_CUBE_MAP_ARRAY 0x900A
#define GL_PROXY_TEXTURE_CUBE_MAP_ARRAY   0x900B
#define GL_SAMPLER_CUBE_MAP_ARRAY         0x900C
#define GL_SAMPLER_CUBE_MAP_ARRAY_SHADOW  0x900D
#define GL_INT_SAMPLER_CUBE_MAP_ARRAY     0x900E
#define GL_UNSIGNED_INT_SAMPLER_CUBE_MAP_ARRAY 0x900F
#define GL_NUM_SHADING_LANGUAGE_VERSIONS  0x82E9
#define GL_VERTEX_ATTRIB_ARRAY_LONG       0x874E
#define GL_DEPTH_COMPONENT32F             0x8CAC
#define GL_DEPTH32F_STENCIL8              0x8CAD
#define GL_FLOAT_32_UNSIGNED_INT_24_8_REV 0x8DAD
#define GL_INVALID_FRAMEBUFFER_OPERATION  0x0506
#define GL_FRAMEBUFFER_ATTACHMENT_COLOR_ENCODING 0x8210
#define GL_FRAMEBUFFER_ATTACHMENT_COMPONENT_TYPE 0x8211
#define GL_FRAMEBUFFER_ATTACHMENT_RED_SIZE 0x8212
#define GL_FRAMEBUFFER_ATTACHMENT_GREEN_SIZE 0x8213
#define GL_FRAMEBUFFER_ATTACHMENT_BLUE_SIZE 0x8214
#define GL_FRAMEBUFFER_ATTACHMENT_ALPHA_SIZE 0x8215
#define GL_FRAMEBUFFER_ATTACHMENT_DEPTH_SIZE 0x8216
#define GL_FRAMEBUFFER_ATTACHMENT_STENCIL_SIZE 0x8217
#define GL_FRAMEBUFFER_DEFAULT            0x8218
#define GL_FRAMEBUFFER_UNDEFINED          0x8219
#define GL_DEPTH_STENCIL_ATTACHMENT       0x821A
#define GL_MAX_RENDERBUFFER_SIZE          0x84E8
#define GL_DEPTH_STENCIL                  0x84F9
#define GL_UNSIGNED_INT_24_8              0x84FA
#define GL_DEPTH24_STENCIL8               0x88F0
#define GL_TEXTURE_STENCIL_SIZE           0x88F1
#define GL_TEXTURE_RED_TYPE               0x8C10
#define GL_TEXTURE_GREEN_TYPE             0x8C11
#define GL_TEXTURE_BLUE_TYPE              0x8C12
#define GL_TEXTURE_ALPHA_TYPE             0x8C13
#define GL_TEXTURE_DEPTH_TYPE             0x8C16
#define GL_UNSIGNED_NORMALIZED            0x8C17
#define GL_FRAMEBUFFER_BINDING            0x8CA6
#define GL_DRAW_FRAMEBUFFER_BINDING       GL_FRAMEBUFFER_BINDING
#define GL_RENDERBUFFER_BINDING           0x8CA7
#define GL_READ_FRAMEBUFFER               0x8CA8
#define GL_DRAW_FRAMEBUFFER               0x8CA9
#define GL_READ_FRAMEBUFFER_BINDING       0x8CAA
#define GL_RENDERBUFFER_SAMPLES           0x8CAB
#define GL_FRAMEBUFFER_ATTACHMENT_OBJECT_TYPE 0x8CD0
#define GL_FRAMEBUFFER_ATTACHMENT_OBJECT_NAME 0x8CD1
#define GL_FRAMEBUFFER_ATTACHMENT_TEXTURE_LEVEL 0x8CD2
#define GL_FRAMEBUFFER_ATTACHMENT_TEXTURE_CUBE_MAP_FACE 0x8CD3
#define GL_FRAMEBUFFER_ATTACHMENT_TEXTURE_LAYER 0x8CD4
#define GL_FRAMEBUFFER_COMPLETE           0x8CD5
#define GL_FRAMEBUFFER_INCOMPLETE_ATTACHMENT 0x8CD6
#define GL_FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT 0x8CD7
#define GL_FRAMEBUFFER_INCOMPLETE_DRAW_BUFFER 0x8CDB
#define GL_FRAMEBUFFER_INCOMPLETE_READ_BUFFER 0x8CDC
#define GL_FRAMEBUFFER_UNSUPPORTED        0x8CDD
#define GL_MAX_COLOR_ATTACHMENTS          0x8CDF
#define GL_COLOR_ATTACHMENT0              0x8CE0
#define GL_COLOR_ATTACHMENT1              0x8CE1
#define GL_COLOR_ATTACHMENT2              0x8CE2
#define GL_COLOR_ATTACHMENT3              0x8CE3
#define GL_COLOR_ATTACHMENT4              0x8CE4
#define GL_COLOR_ATTACHMENT5              0x8CE5
#define GL_COLOR_ATTACHMENT6              0x8CE6
#define GL_COLOR_ATTACHMENT7              0x8CE7
#define GL_COLOR_ATTACHMENT8              0x8CE8
#define GL_COLOR_ATTACHMENT9              0x8CE9
#define GL_COLOR_ATTACHMENT10             0x8CEA
#define GL_COLOR_ATTACHMENT11             0x8CEB
#define GL_COLOR_ATTACHMENT12             0x8CEC
#define GL_COLOR_ATTACHMENT13             0x8CED
#define GL_COLOR_ATTACHMENT14             0x8CEE
#define GL_COLOR_ATTACHMENT15             0x8CEF
#define GL_DEPTH_ATTACHMENT               0x8D00
#define GL_STENCIL_ATTACHMENT             0x8D20
#define GL_FRAMEBUFFER                    0x8D40
#define GL_RENDERBUFFER                   0x8D41
#define GL_RENDERBUFFER_WIDTH             0x8D42
#define GL_RENDERBUFFER_HEIGHT            0x8D43
#define GL_RENDERBUFFER_INTERNAL_FORMAT   0x8D44
#define GL_STENCIL_INDEX1                 0x8D46
#define GL_STENCIL_INDEX4                 0x8D47
#define GL_STENCIL_INDEX8                 0x8D48
#define GL_STENCIL_INDEX16                0x8D49
#define GL_RENDERBUFFER_RED_SIZE          0x8D50
#define GL_RENDERBUFFER_GREEN_SIZE        0x8D51
#define GL_RENDERBUFFER_BLUE_SIZE         0x8D52
#define GL_RENDERBUFFER_ALPHA_SIZE        0x8D53
#define GL_RENDERBUFFER_DEPTH_SIZE        0x8D54
#define GL_RENDERBUFFER_STENCIL_SIZE      0x8D55
#define GL_FRAMEBUFFER_INCOMPLETE_MULTISAMPLE 0x8D56
#define GL_MAX_SAMPLES                    0x8D57
#define GL_FRAMEBUFFER_SRGB               0x8DB9
#define GL_HALF_FLOAT                     0x140B
#define GL_MAP_READ_BIT                   0x0001
#define GL_MAP_WRITE_BIT                  0x0002
#define GL_MAP_INVALIDATE_RANGE_BIT       0x0004
#define GL_MAP_INVALIDATE_BUFFER_BIT      0x0008
#define GL_MAP_FLUSH_EXPLICIT_BIT         0x0010
#define GL_MAP_UNSYNCHRONIZED_BIT         0x0020
#define GL_COMPRESSED_RED_RGTC1           0x8DBB
#define GL_COMPRESSED_SIGNED_RED_RGTC1    0x8DBC
#define GL_COMPRESSED_RG_RGTC2            0x8DBD
#define GL_COMPRESSED_SIGNED_RG_RGTC2     0x8DBE
#define GL_RG                             0x8227
#define GL_RG_INTEGER                     0x8228
#define GL_R8                             0x8229
#define GL_R16                            0x822A
#define GL_RG8                            0x822B
#define GL_RG16                           0x822C
#define GL_R16F                           0x822D
#define GL_R32F                           0x822E
#define GL_RG16F                          0x822F
#define GL_RG32F                          0x8230
#define GL_R8I                            0x8231
#define GL_R8UI                           0x8232
#define GL_R16I                           0x8233
#define GL_R16UI                          0x8234
#define GL_R32I                           0x8235
#define GL_R32UI                          0x8236
#define GL_RG8I                           0x8237
#define GL_RG8UI                          0x8238
#define GL_RG16I                          0x8239
#define GL_RG16UI                         0x823A
#define GL_RG32I                          0x823B
#define GL_RG32UI                         0x823C
#define GL_VERTEX_ARRAY_BINDING           0x85B5
#define GL_UNIFORM_BUFFER                 0x8A11
#define GL_UNIFORM_BUFFER_BINDING         0x8A28
#define GL_UNIFORM_BUFFER_START           0x8A29
#define GL_UNIFORM_BUFFER_SIZE            0x8A2A
#define GL_MAX_VERTEX_UNIFORM_BLOCKS      0x8A2B
#define GL_MAX_GEOMETRY_UNIFORM_BLOCKS    0x8A2C
#define GL_MAX_FRAGMENT_UNIFORM_BLOCKS    0x8A2D
#define GL_MAX_COMBINED_UNIFORM_BLOCKS    0x8A2E
#define GL_MAX_UNIFORM_BUFFER_BINDINGS    0x8A2F
#define GL_MAX_UNIFORM_BLOCK_SIZE         0x8A30
#define GL_MAX_COMBINED_VERTEX_UNIFORM_COMPONENTS 0x8A31
#define GL_MAX_COMBINED_GEOMETRY_UNIFORM_COMPONENTS 0x8A32
#define GL_MAX_COMBINED_FRAGMENT_UNIFORM_COMPONENTS 0x8A33
#define GL_UNIFORM_BUFFER_OFFSET_ALIGNMENT 0x8A34
#define GL_ACTIVE_UNIFORM_BLOCK_MAX_NAME_LENGTH 0x8A35
#define GL_ACTIVE_UNIFORM_BLOCKS          0x8A36
#define GL_UNIFORM_TYPE                   0x8A37
#define GL_UNIFORM_SIZE                   0x8A38
#define GL_UNIFORM_NAME_LENGTH            0x8A39
#define GL_UNIFORM_BLOCK_INDEX            0x8A3A
#define GL_UNIFORM_OFFSET                 0x8A3B
#define GL_UNIFORM_ARRAY_STRIDE           0x8A3C
#define GL_UNIFORM_MATRIX_STRIDE          0x8A3D
#define GL_UNIFORM_IS_ROW_MAJOR           0x8A3E
#define GL_UNIFORM_BLOCK_BINDING          0x8A3F
#define GL_UNIFORM_BLOCK_DATA_SIZE        0x8A40
#define GL_UNIFORM_BLOCK_NAME_LENGTH      0x8A41
#define GL_UNIFORM_BLOCK_ACTIVE_UNIFORMS  0x8A42
#define GL_UNIFORM_BLOCK_ACTIVE_UNIFORM_INDICES 0x8A43
#define GL_UNIFORM_BLOCK_REFERENCED_BY_VERTEX_SHADER 0x8A44
#define GL_UNIFORM_BLOCK_REFERENCED_BY_GEOMETRY_SHADER 0x8A45
#define GL_UNIFORM_BLOCK_REFERENCED_BY_FRAGMENT_SHADER 0x8A46
#define GL_INVALID_INDEX                  0xFFFFFFFFu
#define GL_COPY_READ_BUFFER_BINDING       0x8F36
#define GL_COPY_READ_BUFFER               GL_COPY_READ_BUFFER_BINDING
#define GL_COPY_WRITE_BUFFER_BINDING      0x8F37
#define GL_COPY_WRITE_BUFFER              GL_COPY_WRITE_BUFFER_BINDING
#define GL_DEPTH_CLAMP                    0x864F
#define GL_QUADS_FOLLOW_PROVOKING_VERTEX_CONVENTION 0x8E4C
#define GL_FIRST_VERTEX_CONVENTION        0x8E4D
#define GL_LAST_VERTEX_CONVENTION         0x8E4E
#define GL_PROVOKING_VERTEX               0x8E4F
#define GL_TEXTURE_CUBE_MAP_SEAMLESS      0x884F
#define GL_MAX_SERVER_WAIT_TIMEOUT        0x9111
#define GL_OBJECT_TYPE                    0x9112
#define GL_SYNC_CONDITION                 0x9113
#define GL_SYNC_STATUS                    0x9114
#define GL_SYNC_FLAGS                     0x9115
#define GL_SYNC_FENCE                     0x9116
#define GL_SYNC_GPU_COMMANDS_COMPLETE     0x9117
#define GL_UNSIGNALED                     0x9118
#define GL_SIGNALED                       0x9119
#define GL_ALREADY_SIGNALED               0x911A
#define GL_TIMEOUT_EXPIRED                0x911B
#define GL_CONDITION_SATISFIED            0x911C
#define GL_WAIT_FAILED                    0x911D
#define GL_SYNC_FLUSH_COMMANDS_BIT        0x00000001
#define GL_SAMPLE_POSITION                0x8E50
#define GL_SAMPLE_MASK                    0x8E51
#define GL_SAMPLE_MASK_VALUE              0x8E52
#define GL_MAX_SAMPLE_MASK_WORDS          0x8E59
#define GL_TEXTURE_2D_MULTISAMPLE         0x9100
#define GL_PROXY_TEXTURE_2D_MULTISAMPLE   0x9101
#define GL_TEXTURE_2D_MULTISAMPLE_ARRAY   0x9102
#define GL_PROXY_TEXTURE_2D_MULTISAMPLE_ARRAY 0x9103
#define GL_TEXTURE_BINDING_2D_MULTISAMPLE 0x9104
#define GL_TEXTURE_BINDING_2D_MULTISAMPLE_ARRAY 0x9105
#define GL_TEXTURE_SAMPLES                0x9106
#define GL_TEXTURE_FIXED_SAMPLE_LOCATIONS 0x9107
#define GL_SAMPLER_2D_MULTISAMPLE         0x9108
#define GL_INT_SAMPLER_2D_MULTISAMPLE     0x9109
#define GL_UNSIGNED_INT_SAMPLER_2D_MULTISAMPLE 0x910A
#define GL_SAMPLER_2D_MULTISAMPLE_ARRAY   0x910B
#define GL_INT_SAMPLER_2D_MULTISAMPLE_ARRAY 0x910C
#define GL_UNSIGNED_INT_SAMPLER_2D_MULTISAMPLE_ARRAY 0x910D
#define GL_MAX_COLOR_TEXTURE_SAMPLES      0x910E
#define GL_MAX_DEPTH_TEXTURE_SAMPLES      0x910F
#define GL_MAX_INTEGER_SAMPLES            0x9110
#define GL_SAMPLE_SHADING_ARB             0x8C36
#define GL_MIN_SAMPLE_SHADING_VALUE_ARB   0x8C37
#define GL_TEXTURE_CUBE_MAP_ARRAY_ARB     0x9009
#define GL_TEXTURE_BINDING_CUBE_MAP_ARRAY_ARB 0x900A
#define GL_PROXY_TEXTURE_CUBE_MAP_ARRAY_ARB 0x900B
#define GL_SAMPLER_CUBE_MAP_ARRAY_ARB     0x900C
#define GL_SAMPLER_CUBE_MAP_ARRAY_SHADOW_ARB 0x900D
#define GL_INT_SAMPLER_CUBE_MAP_ARRAY_ARB 0x900E
#define GL_UNSIGNED_INT_SAMPLER_CUBE_MAP_ARRAY_ARB 0x900F
#define GL_MIN_PROGRAM_TEXTURE_GATHER_OFFSET_ARB 0x8E5E
#define GL_MAX_PROGRAM_TEXTURE_GATHER_OFFSET_ARB 0x8E5F
#define GL_MAX_PROGRAM_TEXTURE_GATHER_COMPONENTS_ARB 0x8F9F
#define GL_SHADER_INCLUDE_ARB             0x8DAE
#define GL_NAMED_STRING_LENGTH_ARB        0x8DE9
#define GL_NAMED_STRING_TYPE_ARB          0x8DEA
#define GL_COMPRESSED_RGBA_BPTC_UNORM_ARB 0x8E8C
#define GL_COMPRESSED_SRGB_ALPHA_BPTC_UNORM_ARB 0x8E8D
#define GL_COMPRESSED_RGB_BPTC_SIGNED_FLOAT_ARB 0x8E8E
#define GL_COMPRESSED_RGB_BPTC_UNSIGNED_FLOAT_ARB 0x8E8F
#define GL_SRC1_COLOR                     0x88F9
#define GL_ONE_MINUS_SRC1_COLOR           0x88FA
#define GL_ONE_MINUS_SRC1_ALPHA           0x88FB
#define GL_MAX_DUAL_SOURCE_DRAW_BUFFERS   0x88FC
#define GL_ANY_SAMPLES_PASSED             0x8C2F
#define GL_SAMPLER_BINDING                0x8919
#define GL_RGB10_A2UI                     0x906F
#define GL_TEXTURE_SWIZZLE_R              0x8E42
#define GL_TEXTURE_SWIZZLE_G              0x8E43
#define GL_TEXTURE_SWIZZLE_B              0x8E44
#define GL_TEXTURE_SWIZZLE_A              0x8E45
#define GL_TEXTURE_SWIZZLE_RGBA           0x8E46
#define GL_TIME_ELAPSED                   0x88BF
#define GL_TIMESTAMP                      0x8E28
#define GL_INT_2_10_10_10_REV             0x8D9F
#define GL_DRAW_INDIRECT_BUFFER           0x8F3F
#define GL_DRAW_INDIRECT_BUFFER_BINDING   0x8F43
#define GL_GEOMETRY_SHADER_INVOCATIONS    0x887F
#define GL_MAX_GEOMETRY_SHADER_INVOCATIONS 0x8E5A
#define GL_MIN_FRAGMENT_INTERPOLATION_OFFSET 0x8E5B
#define GL_MAX_FRAGMENT_INTERPOLATION_OFFSET 0x8E5C
#define GL_FRAGMENT_INTERPOLATION_OFFSET_BITS 0x8E5D
#define GL_DOUBLE_VEC2                    0x8FFC
#define GL_DOUBLE_VEC3                    0x8FFD
#define GL_DOUBLE_VEC4                    0x8FFE
#define GL_DOUBLE_MAT2                    0x8F46
#define GL_DOUBLE_MAT3                    0x8F47
#define GL_DOUBLE_MAT4                    0x8F48
#define GL_DOUBLE_MAT2x3                  0x8F49
#define GL_DOUBLE_MAT2x4                  0x8F4A
#define GL_DOUBLE_MAT3x2                  0x8F4B
#define GL_DOUBLE_MAT3x4                  0x8F4C
#define GL_DOUBLE_MAT4x2                  0x8F4D
#define GL_DOUBLE_MAT4x3                  0x8F4E
#define GL_ACTIVE_SUBROUTINES             0x8DE5
#define GL_ACTIVE_SUBROUTINE_UNIFORMS     0x8DE6
#define GL_ACTIVE_SUBROUTINE_UNIFORM_LOCATIONS 0x8E47
#define GL_ACTIVE_SUBROUTINE_MAX_LENGTH   0x8E48
#define GL_ACTIVE_SUBROUTINE_UNIFORM_MAX_LENGTH 0x8E49
#define GL_MAX_SUBROUTINES                0x8DE7
#define GL_MAX_SUBROUTINE_UNIFORM_LOCATIONS 0x8DE8
#define GL_NUM_COMPATIBLE_SUBROUTINES     0x8E4A
#define GL_COMPATIBLE_SUBROUTINES         0x8E4B
#define GL_PATCHES                        0x000E
#define GL_PATCH_VERTICES                 0x8E72
#define GL_PATCH_DEFAULT_INNER_LEVEL      0x8E73
#define GL_PATCH_DEFAULT_OUTER_LEVEL      0x8E74
#define GL_TESS_CONTROL_OUTPUT_VERTICES   0x8E75
#define GL_TESS_GEN_MODE                  0x8E76
#define GL_TESS_GEN_SPACING               0x8E77
#define GL_TESS_GEN_VERTEX_ORDER          0x8E78
#define GL_TESS_GEN_POINT_MODE            0x8E79
#define GL_ISOLINES                       0x8E7A
#define GL_FRACTIONAL_ODD                 0x8E7B
#define GL_FRACTIONAL_EVEN                0x8E7C
#define GL_MAX_PATCH_VERTICES             0x8E7D
#define GL_MAX_TESS_GEN_LEVEL             0x8E7E
#define GL_MAX_TESS_CONTROL_UNIFORM_COMPONENTS 0x8E7F
#define GL_MAX_TESS_EVALUATION_UNIFORM_COMPONENTS 0x8E80
#define GL_MAX_TESS_CONTROL_TEXTURE_IMAGE_UNITS 0x8E81
#define GL_MAX_TESS_EVALUATION_TEXTURE_IMAGE_UNITS 0x8E82
#define GL_MAX_TESS_CONTROL_OUTPUT_COMPONENTS 0x8E83
#define GL_MAX_TESS_PATCH_COMPONENTS      0x8E84
#define GL_MAX_TESS_CONTROL_TOTAL_OUTPUT_COMPONENTS 0x8E85
#define GL_MAX_TESS_EVALUATION_OUTPUT_COMPONENTS 0x8E86
#define GL_MAX_TESS_CONTROL_UNIFORM_BLOCKS 0x8E89
#define GL_MAX_TESS_EVALUATION_UNIFORM_BLOCKS 0x8E8A
#define GL_MAX_TESS_CONTROL_INPUT_COMPONENTS 0x886C
#define GL_MAX_TESS_EVALUATION_INPUT_COMPONENTS 0x886D
#define GL_MAX_COMBINED_TESS_CONTROL_UNIFORM_COMPONENTS 0x8E1E
#define GL_MAX_COMBINED_TESS_EVALUATION_UNIFORM_COMPONENTS 0x8E1F
#define GL_UNIFORM_BLOCK_REFERENCED_BY_TESS_CONTROL_SHADER 0x84F0
#define GL_UNIFORM_BLOCK_REFERENCED_BY_TESS_EVALUATION_SHADER 0x84F1
#define GL_TESS_EVALUATION_SHADER         0x8E87
#define GL_TESS_CONTROL_SHADER            0x8E88
#define GL_TRANSFORM_FEEDBACK             0x8E22
#define GL_TRANSFORM_FEEDBACK_PAUSED      0x8E23
#define GL_TRANSFORM_FEEDBACK_BUFFER_PAUSED GL_TRANSFORM_FEEDBACK_PAUSED
#define GL_TRANSFORM_FEEDBACK_ACTIVE      0x8E24
#define GL_TRANSFORM_FEEDBACK_BUFFER_ACTIVE GL_TRANSFORM_FEEDBACK_ACTIVE
#define GL_TRANSFORM_FEEDBACK_BINDING     0x8E25
#define GL_MAX_TRANSFORM_FEEDBACK_BUFFERS 0x8E70
#define GL_MAX_VERTEX_STREAMS             0x8E71
#define GL_FIXED                          0x140C
#define GL_IMPLEMENTATION_COLOR_READ_TYPE 0x8B9A
#define GL_IMPLEMENTATION_COLOR_READ_FORMAT 0x8B9B
#define GL_LOW_FLOAT                      0x8DF0
#define GL_MEDIUM_FLOAT                   0x8DF1
#define GL_HIGH_FLOAT                     0x8DF2
#define GL_LOW_INT                        0x8DF3
#define GL_MEDIUM_INT                     0x8DF4
#define GL_HIGH_INT                       0x8DF5
#define GL_SHADER_COMPILER                0x8DFA
#define GL_SHADER_BINARY_FORMATS          0x8DF8
#define GL_NUM_SHADER_BINARY_FORMATS      0x8DF9
#define GL_MAX_VERTEX_UNIFORM_VECTORS     0x8DFB
#define GL_MAX_VARYING_VECTORS            0x8DFC
#define GL_MAX_FRAGMENT_UNIFORM_VECTORS   0x8DFD
#define GL_RGB565                         0x8D62
#define GL_PROGRAM_BINARY_RETRIEVABLE_HINT 0x8257
#define GL_PROGRAM_BINARY_LENGTH          0x8741
#define GL_NUM_PROGRAM_BINARY_FORMATS     0x87FE
#define GL_PROGRAM_BINARY_FORMATS         0x87FF
#define GL_VERTEX_SHADER_BIT              0x00000001
#define GL_FRAGMENT_SHADER_BIT            0x00000002
#define GL_GEOMETRY_SHADER_BIT            0x00000004
#define GL_TESS_CONTROL_SHADER_BIT        0x00000008
#define GL_TESS_EVALUATION_SHADER_BIT     0x00000010
#define GL_ALL_SHADER_BITS                0xFFFFFFFF
#define GL_PROGRAM_SEPARABLE              0x8258
#define GL_ACTIVE_PROGRAM                 0x8259
#define GL_PROGRAM_PIPELINE_BINDING       0x825A
#define GL_MAX_VIEWPORTS                  0x825B
#define GL_VIEWPORT_SUBPIXEL_BITS         0x825C
#define GL_VIEWPORT_BOUNDS_RANGE          0x825D
#define GL_LAYER_PROVOKING_VERTEX         0x825E
#define GL_VIEWPORT_INDEX_PROVOKING_VERTEX 0x825F
#define GL_UNDEFINED_VERTEX               0x8260
#define GL_SYNC_CL_EVENT_ARB              0x8240
#define GL_SYNC_CL_EVENT_COMPLETE_ARB     0x8241
#define GL_DEBUG_OUTPUT_SYNCHRONOUS_ARB   0x8242
#define GL_DEBUG_NEXT_LOGGED_MESSAGE_LENGTH_ARB 0x8243
#define GL_DEBUG_CALLBACK_FUNCTION_ARB    0x8244
#define GL_DEBUG_CALLBACK_USER_PARAM_ARB  0x8245
#define GL_DEBUG_SOURCE_API_ARB           0x8246
#define GL_DEBUG_SOURCE_WINDOW_SYSTEM_ARB 0x8247
#define GL_DEBUG_SOURCE_SHADER_COMPILER_ARB 0x8248
#define GL_DEBUG_SOURCE_THIRD_PARTY_ARB   0x8249
#define GL_DEBUG_SOURCE_APPLICATION_ARB   0x824A
#define GL_DEBUG_SOURCE_OTHER_ARB         0x824B
#define GL_DEBUG_TYPE_ERROR_ARB           0x824C
#define GL_DEBUG_TYPE_DEPRECATED_BEHAVIOR_ARB 0x824D
#define GL_DEBUG_TYPE_UNDEFINED_BEHAVIOR_ARB 0x824E
#define GL_DEBUG_TYPE_PORTABILITY_ARB     0x824F
#define GL_DEBUG_TYPE_PERFORMANCE_ARB     0x8250
#define GL_DEBUG_TYPE_OTHER_ARB           0x8251
#define GL_MAX_DEBUG_MESSAGE_LENGTH_ARB   0x9143
#define GL_MAX_DEBUG_LOGGED_MESSAGES_ARB  0x9144
#define GL_DEBUG_LOGGED_MESSAGES_ARB      0x9145
#define GL_DEBUG_SEVERITY_HIGH_ARB        0x9146
#define GL_DEBUG_SEVERITY_MEDIUM_ARB      0x9147
#define GL_DEBUG_SEVERITY_LOW_ARB         0x9148
#define GL_CONTEXT_FLAG_ROBUST_ACCESS_BIT_ARB 0x00000004
#define GL_LOSE_CONTEXT_ON_RESET_ARB      0x8252
#define GL_GUILTY_CONTEXT_RESET_ARB       0x8253
#define GL_INNOCENT_CONTEXT_RESET_ARB     0x8254
#define GL_UNKNOWN_CONTEXT_RESET_ARB      0x8255
#define GL_RESET_NOTIFICATION_STRATEGY_ARB 0x8256
#define GL_NO_RESET_NOTIFICATION_ARB      0x8261
#define GL_UNPACK_COMPRESSED_BLOCK_WIDTH  0x9127
#define GL_UNPACK_COMPRESSED_BLOCK_HEIGHT 0x9128
#define GL_UNPACK_COMPRESSED_BLOCK_DEPTH  0x9129
#define GL_UNPACK_COMPRESSED_BLOCK_SIZE   0x912A
#define GL_PACK_COMPRESSED_BLOCK_WIDTH    0x912B
#define GL_PACK_COMPRESSED_BLOCK_HEIGHT   0x912C
#define GL_PACK_COMPRESSED_BLOCK_DEPTH    0x912D
#define GL_PACK_COMPRESSED_BLOCK_SIZE     0x912E
#define GL_NUM_SAMPLE_COUNTS              0x9380
#define GL_MIN_MAP_BUFFER_ALIGNMENT       0x90BC
#define GL_ATOMIC_COUNTER_BUFFER          0x92C0
#define GL_ATOMIC_COUNTER_BUFFER_BINDING  0x92C1
#define GL_ATOMIC_COUNTER_BUFFER_START    0x92C2
#define GL_ATOMIC_COUNTER_BUFFER_SIZE     0x92C3
#define GL_ATOMIC_COUNTER_BUFFER_DATA_SIZE 0x92C4
#define GL_ATOMIC_COUNTER_BUFFER_ACTIVE_ATOMIC_COUNTERS 0x92C5
#define GL_ATOMIC_COUNTER_BUFFER_ACTIVE_ATOMIC_COUNTER_INDICES 0x92C6
#define GL_ATOMIC_COUNTER_BUFFER_REFERENCED_BY_VERTEX_SHADER 0x92C7
#define GL_ATOMIC_COUNTER_BUFFER_REFERENCED_BY_TESS_CONTROL_SHADER 0x92C8
#define GL_ATOMIC_COUNTER_BUFFER_REFERENCED_BY_TESS_EVALUATION_SHADER 0x92C9
#define GL_ATOMIC_COUNTER_BUFFER_REFERENCED_BY_GEOMETRY_SHADER 0x92CA
#define GL_ATOMIC_COUNTER_BUFFER_REFERENCED_BY_FRAGMENT_SHADER 0x92CB
#define GL_MAX_VERTEX_ATOMIC_COUNTER_BUFFERS 0x92CC
#define GL_MAX_TESS_CONTROL_ATOMIC_COUNTER_BUFFERS 0x92CD
#define GL_MAX_TESS_EVALUATION_ATOMIC_COUNTER_BUFFERS 0x92CE
#define GL_MAX_GEOMETRY_ATOMIC_COUNTER_BUFFERS 0x92CF
#define GL_MAX_FRAGMENT_ATOMIC_COUNTER_BUFFERS 0x92D0
#define GL_MAX_COMBINED_ATOMIC_COUNTER_BUFFERS 0x92D1
#define GL_MAX_VERTEX_ATOMIC_COUNTERS     0x92D2
#define GL_MAX_TESS_CONTROL_ATOMIC_COUNTERS 0x92D3
#define GL_MAX_TESS_EVALUATION_ATOMIC_COUNTERS 0x92D4
#define GL_MAX_GEOMETRY_ATOMIC_COUNTERS   0x92D5
#define GL_MAX_FRAGMENT_ATOMIC_COUNTERS   0x92D6
#define GL_MAX_COMBINED_ATOMIC_COUNTERS   0x92D7
#define GL_MAX_ATOMIC_COUNTER_BUFFER_SIZE 0x92D8
#define GL_MAX_ATOMIC_COUNTER_BUFFER_BINDINGS 0x92DC
#define GL_ACTIVE_ATOMIC_COUNTER_BUFFERS  0x92D9
#define GL_UNIFORM_ATOMIC_COUNTER_BUFFER_INDEX 0x92DA
#define GL_UNSIGNED_INT_ATOMIC_COUNTER    0x92DB
#define GL_VERTEX_ATTRIB_ARRAY_BARRIER_BIT 0x00000001
#define GL_ELEMENT_ARRAY_BARRIER_BIT      0x00000002
#define GL_UNIFORM_BARRIER_BIT            0x00000004
#define GL_TEXTURE_FETCH_BARRIER_BIT      0x00000008
#define GL_SHADER_IMAGE_ACCESS_BARRIER_BIT 0x00000020
#define GL_COMMAND_BARRIER_BIT            0x00000040
#define GL_PIXEL_BUFFER_BARRIER_BIT       0x00000080
#define GL_TEXTURE_UPDATE_BARRIER_BIT     0x00000100
#define GL_BUFFER_UPDATE_BARRIER_BIT      0x00000200
#define GL_FRAMEBUFFER_BARRIER_BIT        0x00000400
#define GL_TRANSFORM_FEEDBACK_BARRIER_BIT 0x00000800
#define GL_ATOMIC_COUNTER_BARRIER_BIT     0x00001000
#define GL_ALL_BARRIER_BITS               0xFFFFFFFF
#define GL_MAX_IMAGE_UNITS                0x8F38
#define GL_MAX_COMBINED_IMAGE_UNITS_AND_FRAGMENT_OUTPUTS 0x8F39
#define GL_IMAGE_BINDING_NAME             0x8F3A
#define GL_IMAGE_BINDING_LEVEL            0x8F3B
#define GL_IMAGE_BINDING_LAYERED          0x8F3C
#define GL_IMAGE_BINDING_LAYER            0x8F3D
#define GL_IMAGE_BINDING_ACCESS           0x8F3E
#define GL_IMAGE_1D                       0x904C
#define GL_IMAGE_2D                       0x904D
#define GL_IMAGE_3D                       0x904E
#define GL_IMAGE_2D_RECT                  0x904F
#define GL_IMAGE_CUBE                     0x9050
#define GL_IMAGE_BUFFER                   0x9051
#define GL_IMAGE_1D_ARRAY                 0x9052
#define GL_IMAGE_2D_ARRAY                 0x9053
#define GL_IMAGE_CUBE_MAP_ARRAY           0x9054
#define GL_IMAGE_2D_MULTISAMPLE           0x9055
#define GL_IMAGE_2D_MULTISAMPLE_ARRAY     0x9056
#define GL_INT_IMAGE_1D                   0x9057
#define GL_INT_IMAGE_2D                   0x9058
#define GL_INT_IMAGE_3D                   0x9059
#define GL_INT_IMAGE_2D_RECT              0x905A
#define GL_INT_IMAGE_CUBE                 0x905B
#define GL_INT_IMAGE_BUFFER               0x905C
#define GL_INT_IMAGE_1D_ARRAY             0x905D
#define GL_INT_IMAGE_2D_ARRAY             0x905E
#define GL_INT_IMAGE_CUBE_MAP_ARRAY       0x905F
#define GL_INT_IMAGE_2D_MULTISAMPLE       0x9060
#define GL_INT_IMAGE_2D_MULTISAMPLE_ARRAY 0x9061
#define GL_UNSIGNED_INT_IMAGE_1D          0x9062
#define GL_UNSIGNED_INT_IMAGE_2D          0x9063
#define GL_UNSIGNED_INT_IMAGE_3D          0x9064
#define GL_UNSIGNED_INT_IMAGE_2D_RECT     0x9065
#define GL_UNSIGNED_INT_IMAGE_CUBE        0x9066
#define GL_UNSIGNED_INT_IMAGE_BUFFER      0x9067
#define GL_UNSIGNED_INT_IMAGE_1D_ARRAY    0x9068
#define GL_UNSIGNED_INT_IMAGE_2D_ARRAY    0x9069
#define GL_UNSIGNED_INT_IMAGE_CUBE_MAP_ARRAY 0x906A
#define GL_UNSIGNED_INT_IMAGE_2D_MULTISAMPLE 0x906B
#define GL_UNSIGNED_INT_IMAGE_2D_MULTISAMPLE_ARRAY 0x906C
#define GL_MAX_IMAGE_SAMPLES              0x906D
#define GL_IMAGE_BINDING_FORMAT           0x906E
#define GL_IMAGE_FORMAT_COMPATIBILITY_TYPE 0x90C7
#define GL_IMAGE_FORMAT_COMPATIBILITY_BY_SIZE 0x90C8
#define GL_IMAGE_FORMAT_COMPATIBILITY_BY_CLASS 0x90C9
#define GL_MAX_VERTEX_IMAGE_UNIFORMS      0x90CA
#define GL_MAX_TESS_CONTROL_IMAGE_UNIFORMS 0x90CB
#define GL_MAX_TESS_EVALUATION_IMAGE_UNIFORMS 0x90CC
#define GL_MAX_GEOMETRY_IMAGE_UNIFORMS    0x90CD
#define GL_MAX_FRAGMENT_IMAGE_UNIFORMS    0x90CE
#define GL_MAX_COMBINED_IMAGE_UNIFORMS    0x90CF
#define GL_TEXTURE_IMMUTABLE_FORMAT       0x912F
#define GL_COMPRESSED_RGBA_ASTC_4x4_KHR   0x93B0
#define GL_COMPRESSED_RGBA_ASTC_5x4_KHR   0x93B1
#define GL_COMPRESSED_RGBA_ASTC_5x5_KHR   0x93B2
#define GL_COMPRESSED_RGBA_ASTC_6x5_KHR   0x93B3
#define GL_COMPRESSED_RGBA_ASTC_6x6_KHR   0x93B4
#define GL_COMPRESSED_RGBA_ASTC_8x5_KHR   0x93B5
#define GL_COMPRESSED_RGBA_ASTC_8x6_KHR   0x93B6
#define GL_COMPRESSED_RGBA_ASTC_8x8_KHR   0x93B7
#define GL_COMPRESSED_RGBA_ASTC_10x5_KHR  0x93B8
#define GL_COMPRESSED_RGBA_ASTC_10x6_KHR  0x93B9
#define GL_COMPRESSED_RGBA_ASTC_10x8_KHR  0x93BA
#define GL_COMPRESSED_RGBA_ASTC_10x10_KHR 0x93BB
#define GL_COMPRESSED_RGBA_ASTC_12x10_KHR 0x93BC
#define GL_COMPRESSED_RGBA_ASTC_12x12_KHR 0x93BD
#define GL_COMPRESSED_SRGB8_ALPHA8_ASTC_4x4_KHR 0x93D0
#define GL_COMPRESSED_SRGB8_ALPHA8_ASTC_5x4_KHR 0x93D1
#define GL_COMPRESSED_SRGB8_ALPHA8_ASTC_5x5_KHR 0x93D2
#define GL_COMPRESSED_SRGB8_ALPHA8_ASTC_6x5_KHR 0x93D3
#define GL_COMPRESSED_SRGB8_ALPHA8_ASTC_6x6_KHR 0x93D4
#define GL_COMPRESSED_SRGB8_ALPHA8_ASTC_8x5_KHR 0x93D5
#define GL_COMPRESSED_SRGB8_ALPHA8_ASTC_8x6_KHR 0x93D6
#define GL_COMPRESSED_SRGB8_ALPHA8_ASTC_8x8_KHR 0x93D7
#define GL_COMPRESSED_SRGB8_ALPHA8_ASTC_10x5_KHR 0x93D8
#define GL_COMPRESSED_SRGB8_ALPHA8_ASTC_10x6_KHR 0x93D9
#define GL_COMPRESSED_SRGB8_ALPHA8_ASTC_10x8_KHR 0x93DA
#define GL_COMPRESSED_SRGB8_ALPHA8_ASTC_10x10_KHR 0x93DB
#define GL_COMPRESSED_SRGB8_ALPHA8_ASTC_12x10_KHR 0x93DC
#define GL_COMPRESSED_SRGB8_ALPHA8_ASTC_12x12_KHR 0x93DD
#define GL_DEBUG_OUTPUT_SYNCHRONOUS       0x8242
#define GL_DEBUG_NEXT_LOGGED_MESSAGE_LENGTH 0x8243
#define GL_DEBUG_CALLBACK_FUNCTION        0x8244
#define GL_DEBUG_CALLBACK_USER_PARAM      0x8245
#define GL_DEBUG_SOURCE_API               0x8246
#define GL_DEBUG_SOURCE_WINDOW_SYSTEM     0x8247
#define GL_DEBUG_SOURCE_SHADER_COMPILER   0x8248
#define GL_DEBUG_SOURCE_THIRD_PARTY       0x8249
#define GL_DEBUG_SOURCE_APPLICATION       0x824A
#define GL_DEBUG_SOURCE_OTHER             0x824B
#define GL_DEBUG_TYPE_ERROR               0x824C
#define GL_DEBUG_TYPE_DEPRECATED_BEHAVIOR 0x824D
#define GL_DEBUG_TYPE_UNDEFINED_BEHAVIOR  0x824E
#define GL_DEBUG_TYPE_PORTABILITY         0x824F
#define GL_DEBUG_TYPE_PERFORMANCE         0x8250
#define GL_DEBUG_TYPE_OTHER               0x8251
#define GL_DEBUG_TYPE_MARKER              0x8268
#define GL_DEBUG_TYPE_PUSH_GROUP          0x8269
#define GL_DEBUG_TYPE_POP_GROUP           0x826A
#define GL_DEBUG_SEVERITY_NOTIFICATION    0x826B
#define GL_MAX_DEBUG_GROUP_STACK_DEPTH    0x826C
#define GL_DEBUG_GROUP_STACK_DEPTH        0x826D
#define GL_BUFFER                         0x82E0
#define GL_SHADER                         0x82E1
#define GL_PROGRAM                        0x82E2
#define GL_QUERY                          0x82E3
#define GL_PROGRAM_PIPELINE               0x82E4
#define GL_SAMPLER                        0x82E6
#define GL_DISPLAY_LIST                   0x82E7
#define GL_MAX_LABEL_LENGTH               0x82E8
#define GL_MAX_DEBUG_MESSAGE_LENGTH       0x9143
#define GL_MAX_DEBUG_LOGGED_MESSAGES      0x9144
#define GL_DEBUG_LOGGED_MESSAGES          0x9145
#define GL_DEBUG_SEVERITY_HIGH            0x9146
#define GL_DEBUG_SEVERITY_MEDIUM          0x9147
#define GL_DEBUG_SEVERITY_LOW             0x9148
#define GL_DEBUG_OUTPUT                   0x92E0
#define GL_CONTEXT_FLAG_DEBUG_BIT         0x00000002
#define GL_COMPUTE_SHADER                 0x91B9
#define GL_MAX_COMPUTE_UNIFORM_BLOCKS     0x91BB
#define GL_MAX_COMPUTE_TEXTURE_IMAGE_UNITS 0x91BC
#define GL_MAX_COMPUTE_IMAGE_UNIFORMS     0x91BD
#define GL_MAX_COMPUTE_SHARED_MEMORY_SIZE 0x8262
#define GL_MAX_COMPUTE_UNIFORM_COMPONENTS 0x8263
#define GL_MAX_COMPUTE_ATOMIC_COUNTER_BUFFERS 0x8264
#define GL_MAX_COMPUTE_ATOMIC_COUNTERS    0x8265
#define GL_MAX_COMBINED_COMPUTE_UNIFORM_COMPONENTS 0x8266
#define GL_MAX_COMPUTE_LOCAL_INVOCATIONS  0x90EB
#define GL_MAX_COMPUTE_WORK_GROUP_COUNT   0x91BE
#define GL_MAX_COMPUTE_WORK_GROUP_SIZE    0x91BF
#define GL_COMPUTE_LOCAL_WORK_SIZE        0x8267
#define GL_UNIFORM_BLOCK_REFERENCED_BY_COMPUTE_SHADER 0x90EC
#define GL_ATOMIC_COUNTER_BUFFER_REFERENCED_BY_COMPUTE_SHADER 0x90ED
#define GL_DISPATCH_INDIRECT_BUFFER       0x90EE
#define GL_DISPATCH_INDIRECT_BUFFER_BINDING 0x90EF
#define GL_COMPUTE_SHADER_BIT             0x00000020
#define GL_TEXTURE_VIEW_MIN_LEVEL         0x82DB
#define GL_TEXTURE_VIEW_NUM_LEVELS        0x82DC
#define GL_TEXTURE_VIEW_MIN_LAYER         0x82DD
#define GL_TEXTURE_VIEW_NUM_LAYERS        0x82DE
#define GL_TEXTURE_IMMUTABLE_LEVELS       0x82DF
#define GL_VERTEX_ATTRIB_BINDING          0x82D4
#define GL_VERTEX_ATTRIB_RELATIVE_OFFSET  0x82D5
#define GL_VERTEX_BINDING_DIVISOR         0x82D6
#define GL_VERTEX_BINDING_OFFSET          0x82D7
#define GL_VERTEX_BINDING_STRIDE          0x82D8
#define GL_MAX_VERTEX_ATTRIB_RELATIVE_OFFSET 0x82D9
#define GL_MAX_VERTEX_ATTRIB_BINDINGS     0x82DA
#define GL_COMPRESSED_RGB8_ETC2           0x9274
#define GL_COMPRESSED_SRGB8_ETC2          0x9275
#define GL_COMPRESSED_RGB8_PUNCHTHROUGH_ALPHA1_ETC2 0x9276
#define GL_COMPRESSED_SRGB8_PUNCHTHROUGH_ALPHA1_ETC2 0x9277
#define GL_COMPRESSED_RGBA8_ETC2_EAC      0x9278
#define GL_COMPRESSED_SRGB8_ALPHA8_ETC2_EAC 0x9279
#define GL_COMPRESSED_R11_EAC             0x9270
#define GL_COMPRESSED_SIGNED_R11_EAC      0x9271
#define GL_COMPRESSED_RG11_EAC            0x9272
#define GL_COMPRESSED_SIGNED_RG11_EAC     0x9273
#define GL_PRIMITIVE_RESTART_FIXED_INDEX  0x8D69
#define GL_ANY_SAMPLES_PASSED_CONSERVATIVE 0x8D6A
#define GL_MAX_ELEMENT_INDEX              0x8D6B
#define GL_MAX_UNIFORM_LOCATIONS          0x826E
#define GL_FRAMEBUFFER_DEFAULT_WIDTH      0x9310
#define GL_FRAMEBUFFER_DEFAULT_HEIGHT     0x9311
#define GL_FRAMEBUFFER_DEFAULT_LAYERS     0x9312
#define GL_FRAMEBUFFER_DEFAULT_SAMPLES    0x9313
#define GL_FRAMEBUFFER_DEFAULT_FIXED_SAMPLE_LOCATIONS 0x9314
#define GL_MAX_FRAMEBUFFER_WIDTH          0x9315
#define GL_MAX_FRAMEBUFFER_HEIGHT         0x9316
#define GL_MAX_FRAMEBUFFER_LAYERS         0x9317
#define GL_MAX_FRAMEBUFFER_SAMPLES        0x9318
#define GL_INTERNALFORMAT_SUPPORTED       0x826F
#define GL_INTERNALFORMAT_PREFERRED       0x8270
#define GL_INTERNALFORMAT_RED_SIZE        0x8271
#define GL_INTERNALFORMAT_GREEN_SIZE      0x8272
#define GL_INTERNALFORMAT_BLUE_SIZE       0x8273
#define GL_INTERNALFORMAT_ALPHA_SIZE      0x8274
#define GL_INTERNALFORMAT_DEPTH_SIZE      0x8275
#define GL_INTERNALFORMAT_STENCIL_SIZE    0x8276
#define GL_INTERNALFORMAT_SHARED_SIZE     0x8277
#define GL_INTERNALFORMAT_RED_TYPE        0x8278
#define GL_INTERNALFORMAT_GREEN_TYPE      0x8279
#define GL_INTERNALFORMAT_BLUE_TYPE       0x827A
#define GL_INTERNALFORMAT_ALPHA_TYPE      0x827B
#define GL_INTERNALFORMAT_DEPTH_TYPE      0x827C
#define GL_INTERNALFORMAT_STENCIL_TYPE    0x827D
#define GL_MAX_WIDTH                      0x827E
#define GL_MAX_HEIGHT                     0x827F
#define GL_MAX_DEPTH                      0x8280
#define GL_MAX_LAYERS                     0x8281
#define GL_MAX_COMBINED_DIMENSIONS        0x8282
#define GL_COLOR_COMPONENTS               0x8283
#define GL_DEPTH_COMPONENTS               0x8284
#define GL_STENCIL_COMPONENTS             0x8285
#define GL_COLOR_RENDERABLE               0x8286
#define GL_DEPTH_RENDERABLE               0x8287
#define GL_STENCIL_RENDERABLE             0x8288
#define GL_FRAMEBUFFER_RENDERABLE         0x8289
#define GL_FRAMEBUFFER_RENDERABLE_LAYERED 0x828A
#define GL_FRAMEBUFFER_BLEND              0x828B
#define GL_READ_PIXELS                    0x828C
#define GL_READ_PIXELS_FORMAT             0x828D
#define GL_READ_PIXELS_TYPE               0x828E
#define GL_TEXTURE_IMAGE_FORMAT           0x828F
#define GL_TEXTURE_IMAGE_TYPE             0x8290
#define GL_GET_TEXTURE_IMAGE_FORMAT       0x8291
#define GL_GET_TEXTURE_IMAGE_TYPE         0x8292
#define GL_MIPMAP                         0x8293
#define GL_MANUAL_GENERATE_MIPMAP         0x8294
#define GL_AUTO_GENERATE_MIPMAP           0x8295
#define GL_COLOR_ENCODING                 0x8296
#define GL_SRGB_READ                      0x8297
#define GL_SRGB_WRITE                     0x8298
#define GL_SRGB_DECODE_ARB                0x8299
#define GL_FILTER                         0x829A
#define GL_VERTEX_TEXTURE                 0x829B
#define GL_TESS_CONTROL_TEXTURE           0x829C
#define GL_TESS_EVALUATION_TEXTURE        0x829D
#define GL_GEOMETRY_TEXTURE               0x829E
#define GL_FRAGMENT_TEXTURE               0x829F
#define GL_COMPUTE_TEXTURE                0x82A0
#define GL_TEXTURE_SHADOW                 0x82A1
#define GL_TEXTURE_GATHER                 0x82A2
#define GL_TEXTURE_GATHER_SHADOW          0x82A3
#define GL_SHADER_IMAGE_LOAD              0x82A4
#define GL_SHADER_IMAGE_STORE             0x82A5
#define GL_SHADER_IMAGE_ATOMIC            0x82A6
#define GL_IMAGE_TEXEL_SIZE               0x82A7
#define GL_IMAGE_COMPATIBILITY_CLASS      0x82A8
#define GL_IMAGE_PIXEL_FORMAT             0x82A9
#define GL_IMAGE_PIXEL_TYPE               0x82AA
#define GL_SIMULTANEOUS_TEXTURE_AND_DEPTH_TEST 0x82AC
#define GL_SIMULTANEOUS_TEXTURE_AND_STENCIL_TEST 0x82AD
#define GL_SIMULTANEOUS_TEXTURE_AND_DEPTH_WRITE 0x82AE
#define GL_SIMULTANEOUS_TEXTURE_AND_STENCIL_WRITE 0x82AF
#define GL_TEXTURE_COMPRESSED_BLOCK_WIDTH 0x82B1
#define GL_TEXTURE_COMPRESSED_BLOCK_HEIGHT 0x82B2
#define GL_TEXTURE_COMPRESSED_BLOCK_SIZE  0x82B3
#define GL_CLEAR_BUFFER                   0x82B4
#define GL_TEXTURE_VIEW                   0x82B5
#define GL_VIEW_COMPATIBILITY_CLASS       0x82B6
#define GL_FULL_SUPPORT                   0x82B7
#define GL_CAVEAT_SUPPORT                 0x82B8
#define GL_IMAGE_CLASS_4_X_32             0x82B9
#define GL_IMAGE_CLASS_2_X_32             0x82BA
#define GL_IMAGE_CLASS_1_X_32             0x82BB
#define GL_IMAGE_CLASS_4_X_16             0x82BC
#define GL_IMAGE_CLASS_2_X_16             0x82BD
#define GL_IMAGE_CLASS_1_X_16             0x82BE
#define GL_IMAGE_CLASS_4_X_8              0x82BF
#define GL_IMAGE_CLASS_2_X_8              0x82C0
#define GL_IMAGE_CLASS_1_X_8              0x82C1
#define GL_IMAGE_CLASS_11_11_10           0x82C2
#define GL_IMAGE_CLASS_10_10_10_2         0x82C3
#define GL_VIEW_CLASS_128_BITS            0x82C4
#define GL_VIEW_CLASS_96_BITS             0x82C5
#define GL_VIEW_CLASS_64_BITS             0x82C6
#define GL_VIEW_CLASS_48_BITS             0x82C7
#define GL_VIEW_CLASS_32_BITS             0x82C8
#define GL_VIEW_CLASS_24_BITS             0x82C9
#define GL_VIEW_CLASS_16_BITS             0x82CA
#define GL_VIEW_CLASS_8_BITS              0x82CB
#define GL_VIEW_CLASS_S3TC_DXT1_RGB       0x82CC
#define GL_VIEW_CLASS_S3TC_DXT1_RGBA      0x82CD
#define GL_VIEW_CLASS_S3TC_DXT3_RGBA      0x82CE
#define GL_VIEW_CLASS_S3TC_DXT5_RGBA      0x82CF
#define GL_VIEW_CLASS_RGTC1_RED           0x82D0
#define GL_VIEW_CLASS_RGTC2_RG            0x82D1
#define GL_VIEW_CLASS_BPTC_UNORM          0x82D2
#define GL_VIEW_CLASS_BPTC_FLOAT          0x82D3
#define GL_UNIFORM                        0x92E1
#define GL_UNIFORM_BLOCK                  0x92E2
#define GL_PROGRAM_INPUT                  0x92E3
#define GL_PROGRAM_OUTPUT                 0x92E4
#define GL_BUFFER_VARIABLE                0x92E5
#define GL_SHADER_STORAGE_BLOCK           0x92E6
#define GL_VERTEX_SUBROUTINE              0x92E8
#define GL_TESS_CONTROL_SUBROUTINE        0x92E9
#define GL_TESS_EVALUATION_SUBROUTINE     0x92EA
#define GL_GEOMETRY_SUBROUTINE            0x92EB
#define GL_FRAGMENT_SUBROUTINE            0x92EC
#define GL_COMPUTE_SUBROUTINE             0x92ED
#define GL_VERTEX_SUBROUTINE_UNIFORM      0x92EE
#define GL_TESS_CONTROL_SUBROUTINE_UNIFORM 0x92EF
#define GL_TESS_EVALUATION_SUBROUTINE_UNIFORM 0x92F0
#define GL_GEOMETRY_SUBROUTINE_UNIFORM    0x92F1
#define GL_FRAGMENT_SUBROUTINE_UNIFORM    0x92F2
#define GL_COMPUTE_SUBROUTINE_UNIFORM     0x92F3
#define GL_TRANSFORM_FEEDBACK_VARYING     0x92F4
#define GL_ACTIVE_RESOURCES               0x92F5
#define GL_MAX_NAME_LENGTH                0x92F6
#define GL_MAX_NUM_ACTIVE_VARIABLES       0x92F7
#define GL_MAX_NUM_COMPATIBLE_SUBROUTINES 0x92F8
#define GL_NAME_LENGTH                    0x92F9
#define GL_TYPE                           0x92FA
#define GL_ARRAY_SIZE                     0x92FB
#define GL_OFFSET                         0x92FC
#define GL_BLOCK_INDEX                    0x92FD
#define GL_ARRAY_STRIDE                   0x92FE
#define GL_MATRIX_STRIDE                  0x92FF
#define GL_IS_ROW_MAJOR                   0x9300
#define GL_ATOMIC_COUNTER_BUFFER_INDEX    0x9301
#define GL_BUFFER_BINDING                 0x9302
#define GL_BUFFER_DATA_SIZE               0x9303
#define GL_NUM_ACTIVE_VARIABLES           0x9304
#define GL_ACTIVE_VARIABLES               0x9305
#define GL_REFERENCED_BY_VERTEX_SHADER    0x9306
#define GL_REFERENCED_BY_TESS_CONTROL_SHADER 0x9307
#define GL_REFERENCED_BY_TESS_EVALUATION_SHADER 0x9308
#define GL_REFERENCED_BY_GEOMETRY_SHADER  0x9309
#define GL_REFERENCED_BY_FRAGMENT_SHADER  0x930A
#define GL_REFERENCED_BY_COMPUTE_SHADER   0x930B
#define GL_TOP_LEVEL_ARRAY_SIZE           0x930C
#define GL_TOP_LEVEL_ARRAY_STRIDE         0x930D
#define GL_LOCATION                       0x930E
#define GL_LOCATION_INDEX                 0x930F
#define GL_IS_PER_PATCH                   0x92E7
#define GL_SHADER_STORAGE_BUFFER          0x90D2
#define GL_SHADER_STORAGE_BUFFER_BINDING  0x90D3
#define GL_SHADER_STORAGE_BUFFER_START    0x90D4
#define GL_SHADER_STORAGE_BUFFER_SIZE     0x90D5
#define GL_MAX_VERTEX_SHADER_STORAGE_BLOCKS 0x90D6
#define GL_MAX_GEOMETRY_SHADER_STORAGE_BLOCKS 0x90D7
#define GL_MAX_TESS_CONTROL_SHADER_STORAGE_BLOCKS 0x90D8
#define GL_MAX_TESS_EVALUATION_SHADER_STORAGE_BLOCKS 0x90D9
#define GL_MAX_FRAGMENT_SHADER_STORAGE_BLOCKS 0x90DA
#define GL_MAX_COMPUTE_SHADER_STORAGE_BLOCKS 0x90DB
#define GL_MAX_COMBINED_SHADER_STORAGE_BLOCKS 0x90DC
#define GL_MAX_SHADER_STORAGE_BUFFER_BINDINGS 0x90DD
#define GL_MAX_SHADER_STORAGE_BLOCK_SIZE  0x90DE
#define GL_SHADER_STORAGE_BUFFER_OFFSET_ALIGNMENT 0x90DF
#define GL_SHADER_STORAGE_BARRIER_BIT     0x2000
#define GL_MAX_COMBINED_SHADER_OUTPUT_RESOURCES GL_MAX_COMBINED_IMAGE_UNITS_AND_FRAGMENT_OUTPUTS
#define GL_DEPTH_STENCIL_TEXTURE_MODE     0x90EA
#define GL_TEXTURE_BUFFER_OFFSET          0x919D
#define GL_TEXTURE_BUFFER_SIZE            0x919E
#define GL_TEXTURE_BUFFER_OFFSET_ALIGNMENT 0x919F
#define GL_VERSION_1_0 1
#define GL_VERSION_1_1 1
#define GL_VERSION_1_2 1
#define GL_VERSION_1_3 1
#define GL_VERSION_1_4 1
#define GL_VERSION_1_5 1
#define GL_VERSION_2_0 1
#define GL_VERSION_2_1 1
#define GL_VERSION_3_0 1
#define GL_VERSION_3_1 1
#define GL_VERSION_3_2 1
#define GL_VERSION_3_3 1
#define GL_VERSION_4_0 1
#define GL_VERSION_4_1 1
#define GL_VERSION_4_2 1
#define GL_VERSION_4_3 1
#define GL_ARB_depth_buffer_float 1
#define GL_ARB_framebuffer_object 1
#define GL_ARB_framebuffer_sRGB 1
#define GL_ARB_half_float_vertex 1
#define GL_ARB_map_buffer_range 1
#define GL_ARB_texture_compression_rgtc 1
#define GL_ARB_texture_rg 1
#define GL_ARB_vertex_array_object 1
#define GL_ARB_uniform_buffer_object 1
#define GL_ARB_copy_buffer 1
#define GL_ARB_depth_clamp 1
#define GL_ARB_draw_elements_base_vertex 1
#define GL_ARB_fragment_coord_conventions 1
#define GL_ARB_provoking_vertex 1
#define GL_ARB_seamless_cube_map 1
#define GL_ARB_sync 1
#define GL_ARB_texture_multisample 1
#define GL_ARB_vertex_array_bgra 1
#define GL_ARB_draw_buffers_blend 1
#define GL_ARB_sample_shading 1
#define GL_ARB_texture_cube_map_array 1
#define GL_ARB_texture_gather 1
#define GL_ARB_texture_query_lod 1
#define GL_ARB_shading_language_include 1
#define GL_ARB_texture_compression_bptc 1
#define GL_ARB_blend_func_extended 1
#define GL_ARB_explicit_attrib_location 1
#define GL_ARB_occlusion_query2 1
#define GL_ARB_sampler_objects 1
#define GL_ARB_shader_bit_encoding 1
#define GL_ARB_texture_rgb10_a2ui 1
#define GL_ARB_texture_swizzle 1
#define GL_ARB_timer_query 1
#define GL_ARB_vertex_type_2_10_10_10_rev 1
#define GL_ARB_draw_indirect 1
#define GL_ARB_gpu_shader5 1
#define GL_ARB_gpu_shader_fp64 1
#define GL_ARB_shader_subroutine 1
#define GL_ARB_tessellation_shader 1
#define GL_ARB_texture_buffer_object_rgb32 1
#define GL_ARB_transform_feedback2 1
#define GL_ARB_transform_feedback3 1
#define GL_ARB_ES2_compatibility 1
#define GL_ARB_get_program_binary 1
#define GL_ARB_separate_shader_objects 1
#define GL_ARB_vertex_attrib_64bit 1
#define GL_ARB_viewport_array 1
#define GL_ARB_cl_event 1
#define GL_ARB_debug_output 1
#define GL_ARB_robustness 1
#define GL_ARB_shader_stencil_export 1
#define GL_ARB_base_instance 1
#define GL_ARB_shading_language_420pack 1
#define GL_ARB_transform_feedback_instanced 1
#define GL_ARB_compressed_texture_pixel_storage 1
#define GL_ARB_conservative_depth 1
#define GL_ARB_internalformat_query 1
#define GL_ARB_map_buffer_alignment 1
#define GL_ARB_shader_atomic_counters 1
#define GL_ARB_shader_image_load_store 1
#define GL_ARB_shading_language_packing 1
#define GL_ARB_texture_storage 1
#define GL_KHR_texture_compression_astc_ldr 1
#define GL_KHR_debug 1
#define GL_ARB_arrays_of_arrays 1
#define GL_ARB_clear_buffer_object 1
#define GL_ARB_compute_shader 1
#define GL_ARB_copy_image 1
#define GL_ARB_texture_view 1
#define GL_ARB_vertex_attrib_binding 1
#define GL_ARB_robustness_isolation 1
#define GL_ARB_ES3_compatibility 1
#define GL_ARB_explicit_uniform_location 1
#define GL_ARB_fragment_layer_viewport 1
#define GL_ARB_framebuffer_no_attachments 1
#define GL_ARB_internalformat_query2 1
#define GL_ARB_invalidate_subdata 1
#define GL_ARB_multi_draw_indirect 1
#define GL_ARB_program_interface_query 1
#define GL_ARB_robust_buffer_access_behavior 1
#define GL_ARB_shader_image_size 1
#define GL_ARB_shader_storage_buffer_object 1
#define GL_ARB_stencil_texturing 1
#define GL_ARB_texture_buffer_range 1
#define GL_ARB_texture_query_levels 1
#define GL_ARB_texture_storage_multisample 1

#endif // CONSTANTS_H
