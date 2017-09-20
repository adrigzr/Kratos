#include <iostream>
#include <cmath>

// Project includes
#include "includes/define.h"
#include "includes/kratos_flags.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "includes/cfd_variables.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "includes/serializer.h"
#include "utilities/geometry_utilities.h"
#include "utilities/math_utils.h"

namespace Kratos {
    template< unsigned int TDim >
    class Eigen {
    public:
        double static Calculate_MatrixTrace3(
            const boost::numeric::ublas::bounded_matrix< double, TDim, TDim> matrix
        )
        {
            return matrix(0,0) + matrix(1,1) + matrix(2,2);
        }

        double static Calculate_MatrixTrace2(
            const boost::numeric::ublas::bounded_matrix< double, TDim, TDim> matrix
        )
        {
            return matrix(0,0) + matrix(1,1);
        }

        void static Calculate_EigenValues2(
            const boost::numeric::ublas::bounded_matrix< double, TDim, TDim> symm_matrix,
            array_1d<double, TDim>& eigen_values
        )
        {
            double temp_1 = 0.5*(symm_matrix(0,0)+symm_matrix(1,1));
            double temp_2 = 0.5*std::sqrt(
                                    std::pow(
                                        symm_matrix(0,0) - symm_matrix(1,1),
                                        2
                                    ) + 
                                    4*std::pow(
                                        symm_matrix(0,1),
                                        2
                                    )
                              );
            eigen_values[0] = temp_1 + temp_2;
            eigen_values[1] = temp_1 - temp_2;
        }

        void static Calculate_EigenValues3(
            const boost::numeric::ublas::bounded_matrix< double, TDim, TDim> symm_matrix,
            array_1d<double, TDim>& eigen_values
        )
        {
            //calculate  the trace of symm_matrix
            double q = Calculate_MatrixTrace3(symm_matrix)/3;
            
            boost::numeric::ublas::bounded_matrix< double, TDim, TDim> identity_matrix;
            identity_matrix(0,0) = 1;
            identity_matrix(0,1) = 0;
            identity_matrix(0,2) = 0;
            identity_matrix(1,0) = 0;
            identity_matrix(1,1) = 1;
            identity_matrix(1,2) = 0;
            identity_matrix(2,0) = 0;
            identity_matrix(2,1) = 0;
            identity_matrix(2,2) = 1;


            boost::numeric::ublas::bounded_matrix< double, TDim, TDim> temp;
            boost::numeric::ublas::bounded_matrix< double, TDim, TDim> result;
            temp = symm_matrix - q*identity_matrix;
            noalias(result) = prod(temp,temp);

            double p = std::sqrt(Calculate_MatrixTrace3(result)/6);
            noalias(temp) = (symm_matrix-q*identity_matrix)/p;
            
            double theta = std::acos(MathUtils<double>::Det(temp)/2)/3;

            double t1 = p*2*std::cos(theta)+q;
            double t2 = p*2*std::cos(theta + 2.0*M_PI/3.0)+q;
            double t3 = p*2*std::cos(theta + 4.0*M_PI/3.0)+q;
            double swap;

            if (t2>t1)
            {
                swap = t1;
                t1 = t2;
                t2 = swap;
            }

            if (t3>t2){
                swap = t2;
                t2 = t3;
                t3 = swap;                
            }

            if (t2>t1)
            {
                swap = t1;
                t1 = t2;
                t2 = swap;                
            }

            eigen_values[0] = t1;
            eigen_values[1] = t2;
            eigen_values[2] = t3;
        }

        void static Calculate_EigenVector2(
            const boost::numeric::ublas::bounded_matrix< double, TDim, TDim> symm_matrix,
            const double eigen_value,
            array_1d<double, TDim>& eigen_vector
        )
        {
            double temp = std::sqrt(
                std::pow(symm_matrix(0,1), 2) + 
                std::pow(symm_matrix(0,0) - eigen_value, 2)
            );
            eigen_vector[0] = symm_matrix(0,1)/temp;
            eigen_vector[1] = -(symm_matrix(0,0)-eigen_value)/temp;
        }

        void static Calculate_EigenVector3(
            const boost::numeric::ublas::bounded_matrix< double, TDim, TDim> symm_matrix,
            const double eigen_value,
            array_1d<double, TDim>& eigen_vector
        )
        {
            array_1d<double, TDim> row0;
            array_1d<double, TDim> row1;
            array_1d<double, TDim> row2;

            row0[0] = symm_matrix(0,0)-eigen_value;
            row0[1] = symm_matrix(0,1);
            row0[2] = symm_matrix(0,2);

            row1[0] = symm_matrix(1,0);
            row1[1] = symm_matrix(1,1)-eigen_value;
            row1[2] = symm_matrix(1,2);

            row2[0] = symm_matrix(2,0);
            row2[1] = symm_matrix(2,1);
            row2[2] = symm_matrix(2,2)-eigen_value;

            
            array_1d<double, TDim> r0xr1;
            array_1d<double, TDim> r0xr2;
            array_1d<double, TDim> r1xr2;

            r0xr1 = MathUtils<double>::CrossProduct(row0, row1);
            r0xr2 = MathUtils<double>::CrossProduct(row0, row2);
            r1xr2 = MathUtils<double>::CrossProduct(row1, row2);            

            double d0 = inner_prod(r0xr1, r0xr1);
            double d1 = inner_prod(r0xr2, r0xr2);
            double d2 = inner_prod(r0xr1, r0xr2);

            double dmax = d0;
            int imax = 0;
            if (d1>dmax) {
                dmax = d1;
                imax = 1;
            }
            if (d2>dmax){
                imax = 2;
            }

            if (imax==0)
                eigen_vector = r0xr1/std::sqrt(d0);
            if (imax==1)
                eigen_vector = r0xr2/std::sqrt(d1);
            if (imax==2)
                eigen_vector = r1xr2/std::sqrt(d2);
        }

        void static Calculate_EigenValues(
            const boost::numeric::ublas::bounded_matrix< double, TDim, TDim > symm_matrix,
            array_1d<double, TDim>& eigen_values
        )
        {
            KRATOS_TRY

            if (TDim == 2)
                Calculate_EigenValues2(symm_matrix, eigen_values);
            else if (TDim == 3)
                Calculate_EigenValues3(symm_matrix, eigen_values);
            else
                KRATOS_THROW_ERROR(std::runtime_error,
                    "Eigen decomposition only supports 2D and 3D.","")

            KRATOS_CATCH("")
        }

        void static Calculate_EigenVector(
            const boost::numeric::ublas::bounded_matrix< double, TDim, TDim > symm_matrix,
            const double eigen_value,
            array_1d<double, TDim>& eigen_vector
        )
        {
            KRATOS_TRY

            if (TDim == 2)
                Calculate_EigenVector2(symm_matrix, eigen_value, eigen_vector);
            else if (TDim == 3)
                Calculate_EigenVector3(symm_matrix, eigen_value, eigen_vector);
            else
                KRATOS_THROW_ERROR(std::runtime_error,
                    "Eigen decomposition only supports 2D and 3D.","")

            KRATOS_CATCH("")
        }        

        void static test_eigen22()
        {
            boost::numeric::ublas::bounded_matrix< double, 2, 2 > symm_matrix;
            symm_matrix(0,0) = -8;
            symm_matrix(0,1) = 2;
            symm_matrix(1,0) = 2;
            symm_matrix(1,1) = -5;

            array_1d<double, 2> eigen_values;
            array_1d<array_1d<double, 2>, 2> eigen_vectors;

            Calculate_EigenValues2(symm_matrix, eigen_values);
            std::cout<<eigen_values<<std::endl;

            Calculate_EigenVector2(symm_matrix, eigen_values[0], eigen_vectors[0]);
            Calculate_EigenVector2(symm_matrix, eigen_values[1], eigen_vectors[1]);

            std::cout<<eigen_vectors<<std::endl;            
        }

        void static test_eigen33()
        {
            boost::numeric::ublas::bounded_matrix< double, 3, 3 > symm_matrix;
            symm_matrix(0,0) = -8;
            symm_matrix(0,1) = 2;
            symm_matrix(0,2) = 3;
            symm_matrix(1,0) = 2;
            symm_matrix(1,1) = -5;
            symm_matrix(1,2) = 6;
            symm_matrix(2,0) = 3;
            symm_matrix(2,1) = 6;
            symm_matrix(2,2) = 9;
            
            array_1d<double, 3> eigen_values;
            array_1d<array_1d<double, 3>, 3> eigen_vectors;

            Calculate_EigenValues3(symm_matrix, eigen_values);
            std::cout<<eigen_values<<std::endl;

            Calculate_EigenVector3(symm_matrix, eigen_values[0], eigen_vectors[0]);
            Calculate_EigenVector3(symm_matrix, eigen_values[1], eigen_vectors[1]);
            Calculate_EigenVector3(symm_matrix, eigen_values[2], eigen_vectors[2]);

            std::cout<<eigen_vectors<<std::endl;

        }
    };
}