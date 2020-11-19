#ifndef RECTANGLE_H
#define RECTANGLE_H

namespace shapes {
    int const aaa= 23;
    class Rectangle {
        public:
            int x0, y0, x1, y1;
            Rectangle();
            Rectangle(int x0, int y0, int x1, int y1);
            ~Rectangle();
            int getArea();
            void getSize(int* width, int* height);
            int get_x();
            void move(int dx, int dy);
    };
}
// cdef extern from "Rectangle.h" namespace "shapes":
// cdef extern from "Rectangle.h" namespace "shapes":
//     cdef cppclass Rectangle:
#endif