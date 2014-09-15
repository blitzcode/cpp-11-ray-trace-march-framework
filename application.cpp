
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
#include "mesh.h"

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
            m_renderer->StopRendering();
            m_renderer->StartRendering();
            break;

        case 's': // Screenshot
            m_renderer->SaveToBMP(("Screenshot_" + DateTimeString() + ".bmp").c_str());
            break;

        case ',': // Decrease sample count
        {
            const uint cnt = m_renderer->GetSampleCount();
            if (cnt > 1)
            {
                m_renderer->StopRendering();
                m_renderer->SetSampleCount(cnt / 2);
                m_renderer->StartRendering();
            }
            break;
        }

        case '.': // Increase sample count
            m_renderer->StopRendering();
            m_renderer->SetSampleCount(m_renderer->GetSampleCount() * 2);
            m_renderer->StartRendering();
            break;

        // Change scene
        case '0': InitializeScene(0); break;
        case '1': InitializeScene(1); break;
        case '2': InitializeScene(2); break;
        case '3': InitializeScene(3); break;
        case '4': InitializeScene(4); break;
        case '5': InitializeScene(5); break;
        case '6': InitializeScene(6); break;
        case '7': InitializeScene(7); break;
        case '8': InitializeScene(8); break;
        case '9': InitializeScene(9); break;
    }
}

void Application::ReshapeFunc(int width, int height)
{
    // Trace("Reshape callback %ix%i", width, height);

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
        ("   %2i FPS | Resolution: %ix%i | Samples [,][.]: %2i | [R]estart Renderer\n"
         "          | [S]creenshot | 2x[ESC] Exit | Scene [0-9]: %i"),
        frames_per_second,
        m_wnd_wdh,
        m_wnd_hgt,
        m_renderer->GetSampleCount(),
        m_cur_scene);
    m_font.DrawStringFixed6x12(4, InvY(12), buf, 0xFF000000);
    m_font.DrawStringFixed6x12(3, InvY(11), buf, 0xFF00FF00);

    // Semi-transparent grey bar to enhance contrast
    glEnable(GL_BLEND);
    glBlendFunc(GL_CONSTANT_ALPHA, GL_ONE_MINUS_CONSTANT_ALPHA);
    glBlendColor(0.0f, 0.0f, 0.0f, 0.5f);
    glBegin(GL_QUADS);
        glColor3f(0.0f, 0.0f, 0.0f);
        glVertex2i(0      , InvY(24));
        glVertex2i(InvX(0), InvY(24));
        glVertex2i(InvX(0), InvY(0));
        glVertex2i(0      , InvY(0));
    glEnd();
    glDisable(GL_BLEND);

    DrawSpinningCube(4, InvY(12), 13);

    m_font.Render();

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

void Application::InitializeScene(uint scene_id)
{
    // Setup mesh and camera parameters
    auto mesh = std::unique_ptr<Mesh>(new Mesh());
    Matrix44f cam_mat;
    float fov = 45.0f;
    switch (scene_id)
    {
        case 0:
            mesh->Read("./meshes/torusknot_column_teapot_plane.dat");
            mesh->NormalizeDimensions();
            cam_mat.BuildLookAtMatrix(Vec3f(-1.00001f, 1.0f, 1.0f), Vec3f(0.0f, -0.2f, 0.0f));
            fov = 51.0f;
            break;

        case 1:
        {
            // Cornell box
            mesh->CornellBox();
            mesh->NormalizeDimensions();

            // Rotated small cube
            Mesh cube;
            cube.Read("./meshes/cube.dat");
            cube.NormalizeDimensions();
            Matrix44f scale, roty, rotx, trans;
            scale.Scaling(0.25f);
            rotx.RotationX(45.0f);
            roty.RotationY(45.0f);
            trans.Translation(0.0f, 0.3f, 0.0f);
            cube.Transform(scale * rotx * roty * trans);
            mesh->AddMesh(cube);

            // Camera
            cam_mat.BuildLookAtMatrix(Vec3f(0.0f, 0.0f, -2.0f), Vec3f(0.0f));
            fov = 51.0f;
            break;
        }

        case 2:
            mesh->Read("./meshes/room_table_chair_tv.dat");
            mesh->NormalizeDimensions();
            cam_mat.BuildLookAtMatrix(Vec3f(-0.47f, 0.15f, -0.3f), Vec3f(1.0f, -0.7f, 0.9f));
            fov = 90.0f;
            break;

        case 3:
        {
            mesh->Read("./meshes/table_chair.dat", true);
            mesh->NormalizeDimensions();
            // Ground plane
            const float quad[4][3]
            {
                { -1.2f, -0.219097f,  1.2f },
                {  1.2f, -0.219097f,  1.2f },
                {  1.2f, -0.219097f, -1.2f },
                { -1.2f, -0.219097f, -1.2f }
            };
            mesh->AddQuad(&quad[0][0]);
            cam_mat.BuildLookAtMatrix(Vec3f(1.001f, 1.002f, -1.0f), Vec3f(0.0f, 0.0f, -0.3f));
            fov = 45.0f;
            break;
        }

        case 4:
        {
            mesh->Read("./meshes/head.dat");
            mesh->NormalizeDimensions();
            Matrix44f roty;
            roty.RotationY(30.0f);
            mesh->Transform(roty);
            cam_mat.BuildLookAtMatrix(Vec3f(0.0f, 0.0f, -1.0f), Vec3f(0.0f));
            fov = 75.0f;
            break;
        }

        case 5:
        {
            mesh->Read("./meshes/room_three_windows_two_columns.dat");
            mesh->NormalizeDimensions();

            Mesh extra_obj;
            extra_obj.Read("./meshes/cat.dat");
            extra_obj.NormalizeDimensions();
            Matrix44f scale, roty, trans;
            scale.Scaling(0.25f);
            roty.RotationY(30.0f);
            trans.Translation(-0.065f, -0.1f, 0.05f);
            extra_obj.Transform(scale * roty * trans);
            mesh->AddMesh(extra_obj);

            // Camera
            cam_mat.BuildLookAtMatrix(Vec3f(-0.2f, 0.0f, -0.33f), Vec3f(0.0f, 0.0f, 0.0f));
            fov = 90.0f;
            break;
        }

        case 6:
        {
            mesh->Read("./meshes/water_surface.dat");
            mesh->NormalizeDimensions();

            Mesh extra_obj;
            extra_obj.Read("./meshes/torus_knot.dat");
            extra_obj.NormalizeDimensions();
            Matrix44f scale, trans;
            scale.Scaling(0.25f);
            trans.Translation(-0.0f, 0.2f, 0.0f);
            extra_obj.Transform(scale * trans);
            mesh->AddMesh(extra_obj);

            // Camera
            cam_mat.BuildLookAtMatrix(Vec3f(-1.0f, 2.0f, -1.0f), Vec3f(0.0f, 0.0f, 0.0f));
            fov = 30.0f;
            break;
        }

        case 7:
        {
            mesh->Read("./meshes/griebel.dat");
            mesh->NormalizeDimensions();

            Mesh extra_obj;
            extra_obj.Read("./meshes/teapot.dat");
            extra_obj.NormalizeDimensions();
            Matrix44f scale, roty,  trans;
            scale.Scaling(0.3f);
            roty.RotationY(90.0f);
            trans.Translation(0.0f, 0.1f, 0.0f);
            extra_obj.Transform(scale * roty * trans);
            mesh->AddMesh(extra_obj);

            // Camera
            cam_mat.BuildLookAtMatrix(Vec3f(0.5, 0.5, 0.0f), Vec3f(0.0f, 0.0f, 0.0f));
            fov = 75.0f;
            break;
        }

        case 8:
        {
            mesh->Read("./meshes/killeroo.dat");
            mesh->NormalizeDimensions();
            // Ground plane
            const float quad[4][3]
            {
                { -0.75f, -0.229267f,  0.75f },
                {  0.75f, -0.229267f,  0.75f },
                {  0.75f, -0.229267f, -0.75f },
                { -0.75f, -0.229267f, -0.75f }
            };
            mesh->AddQuad(&quad[0][0]);
            cam_mat.BuildLookAtMatrix(Vec3f(-1.6f, 1.2f, -1.0f), Vec3f(0.0f, 0.0f, -0.1));
            fov = 30.0f;
            break;
        }

        case 9:
        {
            Matrix44f mat_trans, mat_scale, mat_rotx, mat_roty;

            Mesh dwarf;
            dwarf.Read("./meshes/d3d_dwarf.dat");
            dwarf.NormalizeDimensions();
            mat_trans.Translation(0.0f, 0.500100f, 0.0f);
            dwarf.Transform(mat_trans);
            mesh->AddMesh(dwarf);

            Mesh hand;
            hand.Read("./meshes/hand.dat");
            hand.NormalizeDimensions();
            mat_rotx.RotationX(90.0f);
            mat_roty.RotationY(90.0f);
            mat_trans.Translation(0.7f, 0.490801f, 0.0f);
            hand.Transform(mat_rotx * mat_roty * mat_trans);
            mesh->AddMesh(hand);

            Mesh blob;
            blob.Read("./meshes/blob.dat");
            blob.NormalizeDimensions();
            mat_trans.Translation(-0.8f, 0.278176f, 0.0f);
            mat_scale.Scaling(0.6f);
            blob.Transform(mat_scale * mat_trans);
            mesh->AddMesh(blob);

            const float quad[4][3]
            {
                { -1.5f, 0.0f,  1.0f },
                {  1.5f, 0.0f,  1.0f },
                {  1.5f, 0.0f, -1.0f },
                { -1.5f, 0.0f, -1.0f }
            };
            mesh->AddQuad(&quad[0][0]);

            cam_mat.BuildLookAtMatrix(Vec3f(0.0f, 1.5f, -2.0f), Vec3f(0.0f, 0.0f, 0.0f));
            fov = 60.0f;
            break;
        }

        default:
            return;
    }

    // Create scene
    auto scene = std::unique_ptr<Scene>(new Scene(std::move(mesh), fov, cam_mat));

    // Free existing frame buffer / renderer first
    m_renderer.reset(nullptr);

    // Setup renderer with our new scene
    m_renderer = std::unique_ptr<Renderer>(new Renderer(std::move(scene)));
    m_renderer->Resize(m_wnd_wdh, m_wnd_hgt);

    m_cur_scene = scene_id;
}

int Application::Main(int argc, char **argv)
{
    // Setup OpenGL / GLUT
    SetupGLUT(argc, argv);
    Setup2DOpenGL();
    glClearColor(1.0f, 1.0f, 1.0f, 1.0f);

    if (false)
    {
        // Test grid building code on all meshes
        const char meshes[17][64] =
        {
            "./meshes/blob.dat",
            "./meshes/cat.dat",
            "./meshes/cube.dat",
            "./meshes/griebel.dat",
            "./meshes/hand.dat",
            "./meshes/head.dat",
            "./meshes/killeroo.dat",
            "./meshes/table_chair.dat",
            "./meshes/room_three_windows_two_columns.dat",
            "./meshes/torusknot_column_teapot_plane.dat",
            "./meshes/sphere.dat",
            "./meshes/teapot.dat",
            "./meshes/tiger.dat",
            "./meshes/torus_knot.dat",
            "./meshes/water_surface.dat",
            "./meshes/d3d_dwarf.dat",
            "./meshes/meshes/room_table_chair_tv.dat"
        };
        for (uint i=0; i<16; i++)
        {
            auto test_mesh = std::unique_ptr<Mesh>(new Mesh());
            test_mesh->Read(meshes[i]);
            Grid test_grid(std::move(test_mesh), 64);
        }
    }

    InitializeScene(5);

    glutMainLoop();
    return 0;
}

