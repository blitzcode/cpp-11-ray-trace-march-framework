
#ifndef APPLICATION_H
#define APPLICATION_H

#include <memory>

#include "types.h"
#include "lin_alg.h"
#include "font_rendering.h"
#include "renderer.h"

class Application
{
public:
    int Main(int argc, char **argv);

protected:
    // GLUT callbacks
    void IdleFunc();
    void SpecialKeyCallback(int key, int x, int y);
    void KeyCallback(unsigned char key, int x, int y);
    void ReshapeFunc(int width, int height);
    void DisplayFunc();

    // Callback wrappers
    static Application *m_app;
    static void IdleFunc_();
    static void SpecialKeyCallback_(int key, int x, int y);
    static void KeyCallback_(unsigned char key, int x, int y);
    static void ReshapeFunc_(int width, int height);
    static void DisplayFunc_();

    void Setup2DOpenGL();
    void SetupGLUT(int argc, char **argv);
    void DrawSpinningCube(uint x, uint y, uint width);
    void Shutdown();

    // GUI helpers
    int InvX(int x) { return m_wnd_wdh - x; }
    int InvY(int y) { return m_wnd_hgt - y; }

    FontRendering m_font;
    uint m_wnd_wdh = 640;
    uint m_wnd_hgt = 480;
    std::unique_ptr<Renderer> m_renderer;
};

#endif // APPLICATION_H

