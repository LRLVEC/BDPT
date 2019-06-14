#include <cstdio>
#include <GL/_OpenGL.h>
#include <GL/_Window.h>
#include <_Math.h>
#include <_Time.h>
#include <RayTracing/_RayTracing.h>
#include <GL/_Texture.h>
#include <_STL.h>
#include <_BMP.h>
#include <random>

namespace OpenGL
{
	struct BDPT2 :OpenGL
	{
		struct Renderer :Program
		{
			RayTracing::View view;
			Buffer viewBuffer;
			BufferConfig viewArray;
			VertexAttrib position;

			Renderer(SourceManager* _sm)
				:
				Program(_sm, "Renderer", Vector<VertexAttrib*>{&position}),
				view(),
				viewBuffer(&view),
				viewArray(&viewBuffer, ArrayBuffer),
				position(&viewArray, 0, VertexAttrib::two,
					VertexAttrib::Float, false, sizeof(Math::vec2<float>), 0, 0)
			{
				init();
			}
			void clear()
			{
				glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
				glClear(GL_COLOR_BUFFER_BIT);
			}
			virtual void initBufferData()override
			{
			}
			virtual void run()override
			{
				glDrawArrays(GL_TRIANGLE_FAN, 0, 4);
			}
		};
		struct TracerInit :Computers
		{
			struct CounterData :Buffer::Data
			{
				unsigned int num;
				CounterData()
					:
					num(0)
				{
				}
				virtual void* pointer()override
				{
					return &num;
				}
				virtual unsigned int size()override
				{
					return 4;
				}
			};
			struct TrianglePre :Program
			{
				RayTracing::Model* model;
				TrianglePre(SourceManager* _sm, RayTracing::Model* _model)
					:
					Program(_sm, "TrianglePre"),
					model(_model)
				{
					init();
				}
				virtual void initBufferData()override
				{
				}
				virtual void run()override
				{
					glDispatchCompute((model->geometryNum.data.num.triangleNum + 1023) / 1024, 1, 1);
				}
			};
			struct CirclePre :Program
			{
				RayTracing::Model* model;
				CirclePre(SourceManager* _sm, RayTracing::Model* _model)
					:
					Program(_sm, "CirclePre"),
					model(_model)
				{
					init();
				}
				virtual void initBufferData()override
				{
				}
				virtual void run()override
				{
					glDispatchCompute((model->geometryNum.data.num.circleNum + 1023) / 1024, 1, 1);
				}
			};
			struct DecayOriginCalc :Computers
			{
				struct DecayOriginPre :Program
				{
					DecayOriginPre(SourceManager* _sm)
						:
						Program(_sm, "DecayOriginPre")
					{
						init();
					}
					virtual void initBufferData()override
					{
					}
					virtual void run()override
					{
						glDispatchCompute(8, 1, 1);
					}
				};
				struct DecayOrigin :Program
				{
					DecayOrigin(SourceManager* _sm)
						:
						Program(_sm, "DecayOrigin")
					{
						init();
					}
					virtual void initBufferData()override
					{
					}
					virtual void run()override
					{
						glDispatchCompute(1, 1, 1);
					}
				};

				DecayOriginPre decayOriginPre;
				DecayOrigin decayOrigin;

				DecayOriginCalc(SourceManager* _sm)
					:
					decayOriginPre(_sm),
					decayOrigin(_sm)
				{
				}
				virtual void initBufferData()override
				{
				}
				virtual void run()override
				{
					glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
					decayOriginPre.use();
					decayOriginPre.run();
					glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
					decayOrigin.use();
					decayOrigin.run();
				}
			};
			struct CountLightSources :Program
			{
				struct Info
				{
					unsigned int sourcesNum;
					unsigned int blank[3];
				};
				struct InfoData :Buffer::Data
				{
					Info* info;
					InfoData(Info* _info)
						:
						info(_info)
					{
					}
					virtual void* pointer()override
					{
						return info;
					}
					virtual unsigned int size()override
					{
						return sizeof(Info);
					}
				};
				struct LightSourcesData :Buffer::Data
				{
					struct LightSource
					{
						Math::vec2<unsigned int>index;
						float a;
						float blank;
					};
					Info* info;
					LightSourcesData(Info* _info)
						:
						info(_info)
					{
					}
					virtual void* pointer()override
					{
						return nullptr;
					}
					virtual unsigned int size()override
					{
						return info->sourcesNum * sizeof(LightSource);
					}
				};
				struct IndicesInfo
				{
					unsigned int infoIndex;
					unsigned int lightSourcesIndex;
				};

				Info info;
				InfoData infoData;
				LightSourcesData lightSourcesData;
				Buffer infoBuffer;
				Buffer lightSourcesBuffer;
				BufferConfig infoConfig;
				BufferConfig lightSourcesConfig;

				CountLightSources(SourceManager* _sm, IndicesInfo _indices)
					:
					Program(_sm, "CountLightSources"),
					infoData(&info),
					lightSourcesData(&info),
					infoBuffer(&infoData),
					lightSourcesBuffer(&lightSourcesData),
					infoConfig(&infoBuffer, ShaderStorageBuffer, _indices.infoIndex),
					lightSourcesConfig(&lightSourcesBuffer, ShaderStorageBuffer, _indices.lightSourcesIndex)
				{
					init();
				}
				void dataInit(unsigned int maxSources)
				{
					info.sourcesNum = maxSources;
					infoConfig.dataInit();
					lightSourcesConfig.dataInit();
				}
				virtual void initBufferData()override
				{
				}
				virtual void run()
				{
					glDispatchCompute(1, 1, 1);
				}
			};

			RayTracing::Transform* transform;
			CounterData counterData;
			Buffer test;
			BufferConfig testConfig;
			TrianglePre trianglePre;
			CirclePre circlePre;
			DecayOriginCalc decayOriginCalc;
			CountLightSources countLightSources;
			TracerInit(SourceManager* _sm, RayTracing::FrameScale* _frameScale, RayTracing::Model* _model, RayTracing::Transform* _transform)
				:
				transform(_transform),
				test(&counterData),
				testConfig(&test, ShaderStorageBuffer, 10),
				trianglePre(_sm, _model),
				circlePre(_sm, _model),
				decayOriginCalc(_sm),
				countLightSources(_sm, { 11,12 })
			{
				testConfig.dataInit();
			}
			virtual void initBufferData()override
			{
			}
			virtual void run()override
			{
				if (!trianglePre.model->triangles.GPUUpToDate)
				{
					trianglePre.use();
					trianglePre.run();
				}
				if (!circlePre.model->circles.GPUUpToDate)
				{
					circlePre.use();
					circlePre.run();
				}
				if (transform->moved || trianglePre.model->moved)
				{
					decayOriginCalc.run();
					countLightSources.use();
					countLightSources.run();
					countLightSources.lightSourcesConfig.bind();
					CountLightSources::LightSourcesData::LightSource* a = (CountLightSources::LightSourcesData::LightSource*)
						glMapBuffer(GL_SHADER_STORAGE_BUFFER, GL_READ_ONLY);
					glUnmapBuffer(GL_SHADER_STORAGE_BUFFER);
				}
				trianglePre.model->upToDate();
			}
			void clearCounter()
			{
				counterData.num = 0;
				testConfig.refreshData();
			}
		};

		SourceManager sm;
		bool sizeChanged;
		bool paused;
		RayTracing::FrameScale frameScale;
		RayTracing::Transform transform;
		RayTracing::DecayOriginData decayOriginData;
		RayTracing::Model model;
		Buffer frameSizeBuffer;
		Buffer transBuffer;
		Buffer decayOriginBuffer;
		BufferConfig frameSizeUniform;
		BufferConfig transUniform;
		BufferConfig decayOriginStorage;
		BMPData lightSource;
		BMPCubeData cubeData;
		STL stl;
		Texture textures;
		Texture imageTexture;
		TextureCube cube;
		TextureConfig<TextureStorage3D>textureConfig;
		TextureConfig<TextureStorage2D>imageConfig;
		TracerInit tracerInit;
		Renderer renderer;

		BDPT2(Math::vec2<unsigned int> const& _size)
			:
			sm(),
			sizeChanged(true),
			paused(true),
			frameScale(_size),
			transform({ {60.0},{0.002,0.9,0.001},{0.05},{0,0,0},700.0 }),
			model({ {1,2},{3},{4},{5},{6},{7},{3},{9} }),
			frameSizeBuffer(&frameScale),
			transBuffer(&transform.bufferData),
			decayOriginBuffer(&decayOriginData),
			frameSizeUniform(&frameSizeBuffer, UniformBuffer, 0),
			transUniform(&transBuffer, UniformBuffer, 1),
			decayOriginStorage(&decayOriginBuffer, ShaderStorageBuffer, 8),
			lightSource("resources/lightSource.bmp"),
			cubeData("resources/room/"),
			stl(sm.folder.find("resources/Box1.stl").readSTL()),
			textures(&lightSource, 1),
			imageTexture(nullptr, 0),
			cube(&cubeData, 2, RGBA32f, 1, cubeData.bmp[0].header.width, cubeData.bmp[0].header.height),
			textureConfig(&textures, Texture2DArray, RGBA32f, 1, lightSource.bmp.header.width, lightSource.bmp.header.height, 4),
			imageConfig(&imageTexture, Texture2D, RGBA32f, 1, _size.data[0], _size.data[1]),
			tracerInit(&sm, &frameScale, &model, &transform),
			renderer(&sm)
			//frameBuffer(0),
			//renderBuffer(0)
		{
			glBindImageTexture(2, imageTexture.texture, 0, GL_FALSE, 0, GL_READ_WRITE, GL_RGBA32F);

			textureConfig.dataRefresh(0, TextureInputBGRInt, TextureInputUByte, 0, 0, 0, lightSource.bmp.header.width, lightSource.bmp.header.height, 1);
			textureConfig.parameteri(TextureParameter::TextureMagFilter, TextureParameter::MinFilter_Nearest);
			cube.dataInit(0, TextureInputBGRInt, TextureInputUByte);
			renderer.use();
			textures.bindUnit();
			cube.bindUnit();
			imageTexture.bindUnit();
			glEnable(GL_TEXTURE_CUBE_MAP_SEAMLESS);
			/*glDisable(GL_DEPTH_TEST);
			glEnable(GL_BLEND);
			glBlendEquation(GL_FUNC_ADD);
			glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
			glDrawBuffer(GL_FRONT_AND_BACK);*/
			/*
			model.circles.data.circles +=
			{
				{
					{0, 0, 1, -0.95},
					{ 0,0.5,0.95 }, 0.01,
					{ 1,0,0 },
					{
						0,-1,
						0,-1,
						0,-1,
						40,-1,
						0,
						1
					}
				}
			};
			model.spheres.data.spheres +=
			{
				{
					{0.5, -0.3, 0.1, 0.1},
						0, 0,
					{
						1,-1,
						1,-1,
						0,-1,
						0,-1,
						0,
						1.5
					}
				}
			};*/
			/*model.pointLights.data.pointLights +=
			{
				{
					{0.2, 0.2, 0.2},
					{ 0,0,-0.95 }
				}
			};*/
			unsigned int k(model.triangles.trianglesOrigin.trianglesOrigin.length);
			float glow = 100;
			model.addSTL
			(
				stl,
				{
					0,-1,
					0,-1,
					1,-1,
					0,-1,
					0,//{ -0.03,0,-0.03 },
					1
				},
				stl.triangles.length
			);
			for (int c0(0); c0 < stl.triangles.length - k; ++c0)
			{
				::printf("%d: ", c0);
				stl.triangles[c0].print();
			}
			model.triangles.trianglesOrigin.trianglesOrigin.data[0 + k].color.d = 0;
			model.triangles.trianglesOrigin.trianglesOrigin.data[0 + k].color.g = glow;
			model.triangles.trianglesOrigin.trianglesOrigin.data[0 + k].uv1 = { 1,0 };
			model.triangles.trianglesOrigin.trianglesOrigin.data[0 + k].uv2 = { 0,1 };
			model.triangles.trianglesOrigin.trianglesOrigin.data[0 + k].uv3 = { 1,1 };
			//model.triangles.trianglesOrigin.trianglesOrigin.data[30 + k].color.texG = 0;
			model.triangles.trianglesOrigin.trianglesOrigin.data[1 + k].color.d = 0;
			model.triangles.trianglesOrigin.trianglesOrigin.data[1 + k].color.g = glow;
			model.triangles.trianglesOrigin.trianglesOrigin.data[1 + k].uv1 = { 1,1 };
			model.triangles.trianglesOrigin.trianglesOrigin.data[1 + k].uv2 = { 0,1 };
			model.triangles.trianglesOrigin.trianglesOrigin.data[1 + k].uv3 = { 0,0 };
			//model.triangles.trianglesOrigin.trianglesOrigin.data[31 + k].color.texG = 0;
			model.triangles.trianglesOrigin.trianglesOrigin.data[2 + k].color.d = 0.1;
			model.triangles.trianglesOrigin.trianglesOrigin.data[3 + k].color.d = 0.1;
			/*for (int c0(12); c0 < 24; ++c0)
			{
				model.triangles.trianglesOrigin.trianglesOrigin.data[c0 + k].color.d = 0;
				model.triangles.trianglesOrigin.trianglesOrigin.data[c0 + k].color.r = 1;
				model.triangles.trianglesOrigin.trianglesOrigin.data[c0 + k].color.t = 1;
				model.triangles.trianglesOrigin.trianglesOrigin.data[c0 + k].color.n = 1.5;
			}*/
			model.triangles.trianglesOrigin.trianglesOrigin.data[24 + k].color.d = { 1,0,0 };
			model.triangles.trianglesOrigin.trianglesOrigin.data[25 + k].color.d = { 1,0,0 };
			model.triangles.trianglesOrigin.trianglesOrigin.data[28 + k].color.d = { 0,1,0 };
			model.triangles.trianglesOrigin.trianglesOrigin.data[29 + k].color.d = { 0,1,0 };

			model.triangles.numChanged = true;
			model.spheres.numChanged = true;
			model.circles.numChanged = true;
			model.cylinders.numChanged = true;
			model.cones.numChanged = true;
			model.pointLights.numChanged = true;
		}
		virtual void init(FrameScale const& _size) override
		{
			glViewport(0, 0, _size.w, _size.h);
			transform.init(_size);
			frameScale.scale = { (unsigned int)_size.w,(unsigned int)_size.h };
			renderer.viewArray.dataInit();
			frameSizeUniform.dataInit();
			sizeChanged = false;
			transUniform.dataInit();
			decayOriginStorage.dataInit();
			model.dataInit();
			tracerInit.countLightSources.dataInit(4);
			//glBindFramebuffer(GL_FRAMEBUFFER, frameBuffer);
			//glViewport(0, 0, frameScale.scale.data[0], frameScale.scale.data[1]);
			//renderer.clear();
		}
		virtual void run() override
		{
			if (sizeChanged)
			{
				glViewport(0, 0, frameScale.scale.data[0], frameScale.scale.data[1]);
				//glBindRenderbuffer(GL_RENDERBUFFER, renderBuffer);
				//glRenderbufferStorage(GL_RENDERBUFFER, GL_RGBA32F, frameScale.scale.data[0], frameScale.scale.data[1]);
				imageConfig.width = frameScale.scale.data[0];
				imageConfig.height = frameScale.scale.data[1];
				glDeleteTextures(1, &imageTexture.texture);
				imageTexture.create();
				imageConfig.bind();
				imageConfig.allocData();
				glBindImageTexture(2, imageTexture.texture, 0, GL_FALSE, 0, GL_READ_WRITE, GL_RGBA32F);
				//renderer.use();
				//imageTexture.bindUnit();
				tracerInit.clearCounter();
				frameSizeUniform.refreshData();
				transform.bufferData.trans.times = 0;
				sizeChanged = false;
			}
			transform.operate();
			if (transform.updated)
			{
				tracerInit.clearCounter();
				transform.updated = false;
			}
			transUniform.refreshData();
			glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
			tracerInit.run();
			glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
			renderer.use();
			renderer.run();
		}
		virtual void frameSize(int _w, int _h) override
		{
			frameScale.scale = { unsigned int(_w),unsigned int(_h) };
			transform.persp.y = _h;
			transform.persp.updated = true;
			sizeChanged = true;
		}
		virtual void framePos(int, int) override
		{
		}
		virtual void frameFocus(int) override
		{
		}
		virtual void mouseButton(int _button, int _action, int _mods) override
		{
			switch (_button)
			{
				case GLFW_MOUSE_BUTTON_LEFT:transform.mouse.refreshButton(0, _action); break;
				case GLFW_MOUSE_BUTTON_MIDDLE:transform.mouse.refreshButton(1, _action); break;
				case GLFW_MOUSE_BUTTON_RIGHT:transform.mouse.refreshButton(2, _action); break;
			}
		}
		virtual void mousePos(double _x, double _y) override
		{
			transform.mouse.refreshPos(_x, _y);
		}
		virtual void mouseScroll(double _x, double _y) override
		{
			if (_y != 0.0)
				transform.scroll.refresh(_y);
		}
		virtual void key(GLFWwindow* _window, int _key, int _scancode, int _action, int _mods) override
		{
			switch (_key)
			{
				case GLFW_KEY_ESCAPE:
					if (_action == GLFW_PRESS)
						glfwSetWindowShouldClose(_window, true);
					break;
				case GLFW_KEY_A:transform.key.refresh(0, _action); break;
				case GLFW_KEY_D:transform.key.refresh(1, _action); break;
				case GLFW_KEY_W:transform.key.refresh(2, _action); break;
				case GLFW_KEY_S:transform.key.refresh(3, _action); break;
				case GLFW_KEY_P: if (_action == GLFW_PRESS)paused = !paused;
			}
		}
	};
}

int main()
{
	OpenGL::OpenGLInit init(4, 5);
	Window::Window::Data winPara
	{
		"BidirectionalPathTracing",
		{
			{300,300},
			true, false,
		}
	};
	Window::WindowManager wm(winPara);
	OpenGL::BDPT2 bdpt2({ 300,300 });
	wm.init(0, &bdpt2);
	glfwSwapInterval(0);
	FPS fps;
	fps.refresh();
	::printf("FPS:\n");
	while (!wm.close())
	{
		wm.pullEvents();
		wm.render();
		glFinish();
		wm.swapBuffers();
		fps.refresh();
		bdpt2.tracerInit.testConfig.bind();
		::printf("\r%.2lf    %f", fps.fps, bdpt2.transform.bufferData.trans.times);
	}
	return 0;
}
