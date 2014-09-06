
#include "application.h"

#include <cstdlib>
#include <cassert>
#include <chrono>
#include <thread>

#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>

#include "timer.h"
#include "trace.h"

void Application::IdleFunc()
{
    static double last_tick = 0.0;
    const double  max_fps   = 30.0;

    // Frame limiter
    while (true)
    {
        const double cur_tick       = TimerGetTick();
        const double elapsed        = cur_tick - last_tick;
        const double elapsed_needed = 1.0 / max_fps;

        if (elapsed < elapsed_needed)
        {
            std::chrono::microseconds sleep_amount
                (uint((elapsed_needed - elapsed) * 1000000.0 * 0.9));
            std::this_thread::sleep_for(sleep_amount);

            continue;
        }

        last_tick = cur_tick;
        break;
    }

    glutPostRedisplay();
}

void Application::SpecialKeyCallback(int key, int x, int y)
{
}

void Application::KeyCallback(unsigned char key, int x, int y)
{
    switch (key)
    {
        case 27: // Escape, quit application
        {
            static double last_esc_press = -1.0;
            const double cur_tick = TimerGetTick();
            if (cur_tick - last_esc_press < 0.5) // Need to press ESC twice, quickly
                Shutdown();
            last_esc_press = cur_tick;
            break;
        }

        case 'r': // (Re)start rendering
            m_renderer->RestartRendering();
            break;
    }
}

void Application::ReshapeFunc(int width, int height)
{
    Trace("Reshape callback %ix%i", width, height);

    m_wnd_wdh = uint(width);
    m_wnd_hgt = uint(height);

    // Adjust OpenGL
    glViewport(0, 0, width, height);
    Setup2DOpenGL();

    m_renderer->Resize(width, height);
}

void Application::Setup2DOpenGL()
{
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0.0, GLdouble(m_wnd_wdh), 0.0, GLdouble(m_wnd_hgt));
}

void Application::DisplayFunc()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // FPS counter
    static uint frames_per_second = 0;
    {
        static double last_tick = 0;
        static uint num_frames = 0;
        const double cur_tick = TimerGetTick();
        if (cur_tick - last_tick > 1.0)
        {
            last_tick = cur_tick;
            frames_per_second = num_frames;
            num_frames = 0;
        }
        num_frames++;
    }

    m_renderer->Draw(0, 0, m_wnd_wdh, m_wnd_hgt);

    // Frame counter and help text
    char buf[256];
    std::snprintf(
        buf,
        sizeof(buf),
        "%i FPS | 2x[ESC] Exit | [R]estart Renderer",
        frames_per_second);
    m_font.DrawStringFixed6x12(19, InvY(13), buf, 0xFF000000);
    m_font.DrawStringFixed6x12(18, InvY(12), buf, 0xFF00FF00);

    // Semi-transparent grey bar to enhance contrast
    glEnable(GL_BLEND);
    glBlendFunc(GL_CONSTANT_ALPHA, GL_ONE_MINUS_CONSTANT_ALPHA);
    glBlendColor(0.0f, 0.0f, 0.0f, 0.5f);
    glBegin(GL_QUADS);
        glColor3f(0.0f, 0.0f, 0.0f);
        glVertex2i(0      , InvY(14));
        glVertex2i(InvX(0), InvY(14));
        glVertex2i(InvX(0), InvY(0));
        glVertex2i(0      , InvY(0));
    glEnd();
    glDisable(GL_BLEND);

    DrawSpinningCube(2, InvY(13), 13);

    m_font.Render(false);

    // Error checking
#ifndef NDEBUG
    const GLenum error = glGetError();
#endif // NDEBUG
    assert(error != GL_INVALID_ENUM);
    assert(error != GL_INVALID_VALUE);
    assert(error != GL_INVALID_OPERATION);
    assert(error != GL_STACK_OVERFLOW);
    assert(error != GL_STACK_UNDERFLOW);
    assert(error != GL_OUT_OF_MEMORY);
    assert(error != GL_TABLE_TOO_LARGE);
    assert(error == GL_NO_ERROR);

    glutSwapBuffers();
}

void Application::DrawSpinningCube(uint x, uint y, uint width)
{
    // Draw a spinning cube as an indicator that the GUI is still updating

    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    glOrtho(0, m_wnd_wdh, 0, m_wnd_hgt, 1.0f, width * 2);

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    glTranslatef(x + width / 2.0f, y + width / 2.0f, - int(width));
    glScalef(width / 4.0f, width / 4.0f, width / 4.0f);
    const double cur_tick = TimerGetTick();
    glRotatef(360.0f * float(cur_tick / 4.0 - std::floor(cur_tick / 4.0)),
        0.0f, 1.0f, 1.0f);

    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glDisable(GL_CULL_FACE);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_NORMALIZE);

    // TODO: Winding isn't consistent on these
    glBegin(GL_QUADS);
        glNormal3f( 0.0f, -1.0f,  0.0f);
        glVertex3f(-1.0f, -1.0f, -1.0f);
        glVertex3f( 1.0f, -1.0f, -1.0f);
        glVertex3f( 1.0f, -1.0f,  1.0f);
        glVertex3f(-1.0f, -1.0f,  1.0f);
        glNormal3f( 0.0f,  1.0f,  0.0f);
        glVertex3f(-1.0f,  1.0f, -1.0f);
        glVertex3f( 1.0f,  1.0f, -1.0f);
        glVertex3f( 1.0f,  1.0f,  1.0f);
        glVertex3f(-1.0f,  1.0f,  1.0f);
        glNormal3f( 0.0f,  0.0f, -1.0f);
        glVertex3f(-1.0f, -1.0f, -1.0f);
        glVertex3f( 1.0f, -1.0f, -1.0f);
        glVertex3f( 1.0f,  1.0f, -1.0f);
        glVertex3f(-1.0f,  1.0f, -1.0f);
        glNormal3f( 0.0f,  0.0f,  1.0f);
        glVertex3f(-1.0f, -1.0f,  1.0f);
        glVertex3f( 1.0f, -1.0f,  1.0f);
        glVertex3f( 1.0f,  1.0f,  1.0f);
        glVertex3f(-1.0f,  1.0f,  1.0f);
        glNormal3f(-1.0f,  0.0f,  0.0f);
        glVertex3f(-1.0f, -1.0f, -1.0f);
        glVertex3f(-1.0f,  1.0f, -1.0f);
        glVertex3f(-1.0f,  1.0f,  1.0f);
        glVertex3f(-1.0f, -1.0f,  1.0f);
        glNormal3f( 1.0f,  0.0f,  0.0f);
        glVertex3f( 1.0f, -1.0f, -1.0f);
        glVertex3f( 1.0f,  1.0f, -1.0f);
        glVertex3f( 1.0f,  1.0f,  1.0f);
        glVertex3f( 1.0f, -1.0f,  1.0f);
    glEnd();

    glDisable(GL_NORMALIZE);
    glDisable(GL_LIGHTING);
    glDisable(GL_DEPTH_TEST);
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();
}

// Callback wrappers
Application *Application::m_app = nullptr;
void Application::IdleFunc_()
    { m_app->IdleFunc(); }
void Application::SpecialKeyCallback_(int key, int x, int y)
    { m_app->SpecialKeyCallback(key, x, y); }
void Application::KeyCallback_(unsigned char key, int x, int y)
    { m_app->KeyCallback(key, x, y); }
void Application::ReshapeFunc_(int width, int height)
    { m_app->ReshapeFunc(width, height); }
void Application::DisplayFunc_()
    { m_app->DisplayFunc(); }

void Application::SetupGLUT(int argc, char **argv)
{
    // Initialize GLUT
    glutInit(&argc, argv);

    // Window
    glutInitWindowSize(m_wnd_wdh, m_wnd_hgt);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
    glutCreateWindow("Viewer");

    // Callbacks
    m_app = this;
    glutDisplayFunc(DisplayFunc_);
    glutSpecialFunc(SpecialKeyCallback_);
    glutKeyboardFunc(KeyCallback_);
    glutReshapeFunc(ReshapeFunc_);
    glutIdleFunc(IdleFunc_);
}

void Application::Shutdown()
{
    Trace("Shutting down");

    // Since we can never leave GLUTs main loop, no dtors will ever be called. Destroy
    // some objects manually to at least test that their destruction behavior is working
    m_renderer.reset(nullptr);

    exit(0);
}

int Application::Main(int argc, char **argv)
{
    // Setup OpenGL / GLUT
    SetupGLUT(argc, argv);
    Setup2DOpenGL();
    glClearColor(1.0f, 1.0f, 1.0f, 1.0f);

    m_renderer = std::unique_ptr<Renderer>(new Renderer(m_wnd_wdh, m_wnd_hgt));

    glutMainLoop();
    return 0;
}

