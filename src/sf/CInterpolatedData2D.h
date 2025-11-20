#ifndef C_INTERPOLATED_DATA_2D_H
#define C_INTERPOLATED_DATA_2D_H

#include <cmath>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include "GConstants.h"

///   Class handling 2D data tables

///   The input format is the following:
///
///   first the header:
///   xRes xStart xStop xStep
///   yRes yStart yStop yStep
///
///   then xRes blocks of data:
///   the x value, and the current y value and the function value repeated yRes times

class CInterpolatedData2D
{
public:
   CInterpolatedData2D(const std::string inputData,
                       const double i_abscissaUnit,
                       const double i_ordinateUnit,
                       const double i_valueUnit,
                       const bool i_histogram)
       : m_histogram(i_histogram), m_valid(false)
   {
        std::ifstream dataFile( inputData.c_str() );

        if (dataFile.fail())
        {
            std::cerr<<"Indispensable file '"<<inputData<<"' not found"<<std::endl;
            return;
        }

        double dummy(0.0), x(0.0), y(0.0);

        dataFile>>m_xRes>>m_xStart>>m_xStop>>m_xStep;
        dataFile>>m_yRes>>m_yStart>>m_yStop>>m_yStep;

        m_xMinI = m_xStart;
        m_yMinI = m_yStart;

        if ( m_histogram )
        {
            m_xStart -= 0.5*m_xStep;
            m_xStop  += 0.5*m_xStep;

            m_yStart -= 0.5*m_yStep;
            m_yStop  += 0.5*m_yStep;
        }

        m_xStart *= i_abscissaUnit;
        m_xStep  *= i_abscissaUnit;
        m_xStop  *= i_abscissaUnit;

        m_yStart *= i_ordinateUnit;
        m_yStep  *= i_ordinateUnit;
        m_yStop  *= i_ordinateUnit;

        m_dataTable = new double *[m_xRes + 1]; // 1 is added to avoid problems with interpolation at the border
        m_dataTable[m_xRes] = new double[m_yRes + 1]; // 1 is added to avoid problems with interpolation at the border
        m_xDistributionTable = new double [m_xRes];

        double sum(0.0); // For the normalization

        for (int lCnt(0); lCnt < m_xRes; ++lCnt)
        {
            m_dataTable[lCnt] = new double[m_yRes + 1]; // 1 is added to avoid problems with interpolation at the border
            m_xDistributionTable[lCnt] = 0.0;

            dataFile>>x; // Read the current value of x

            for (int kCnt(0); kCnt < m_yRes; ++kCnt)
            {
                dataFile>>y; // Read the current value of y
                dataFile>>dummy; // Read the current function value
                m_dataTable[lCnt][kCnt] = dummy * i_valueUnit;
                m_xDistributionTable[lCnt] += m_dataTable[lCnt][kCnt] * m_yStep;
            }

            sum += x * x * m_xDistributionTable[lCnt];

            m_dataTable[lCnt][m_yRes] = 0.0;
        }

        for (int kCnt(0); kCnt < m_yRes + 1; ++kCnt)
            m_dataTable[m_xRes][kCnt] = 0.0;

        m_normalization = 4 * Pi * sum * m_xStep;

        // Check for zero normalization to prevent division by zero
        if (m_normalization == 0.0)
        {
            std::cerr << "Normalization factor is zero in CInterpolatedData2D for file '" << inputData << "'" << std::endl;
            return;
        }

        m_valid = true; // Successfully initialized

        // std::cout<<"The file '"<<inputData<<"' has been loaded"<<std::endl;
    };

    ~CInterpolatedData2D()
    {
        for (int iCnt(0); iCnt < m_xRes + 1; ++iCnt)
            delete [] m_dataTable[iCnt];

        delete [] m_dataTable;
        delete [] m_xDistributionTable;
    };

    // Check if the object is valid
    bool isValid() const { return m_valid; }

    double withoutI(const double x, const double y) const
    {
        if ( x < m_xStart or y < m_yStart or x > m_xStop or y > m_yStop )
            return 0.0;

        const double xDiv( (x - m_xStart)/m_xStep );
        const double yDiv( (y - m_yStart)/m_yStep );

        const int nx( static_cast<int>(round(xDiv)) );
        const int ny( static_cast<int>(round(yDiv)) );

        const double c00( m_dataTable[ nx ][ ny ] );

        return c00;
    }

    double linearI(const double x, const double y) const
    {
        if ( x < m_xStart or y < m_yStart or x > m_xStop or y > m_yStop )
        {
            std::cout<<x<<" "<<y<<std::endl;
            std::cout<<"You are requesting variables from a range not covered by the table"<<std::endl;
            std::cout<<m_xStart<<" < x < "<<m_xStop<<", "<<m_yStart<<" < y < "<<m_yStop<<std::endl;
            // system("pause");
            return 0.0;
        }

        const double xDiv( (x - m_xMinI)/m_xStep );
        const double yDiv( (y - m_yMinI)/m_yStep );

        const int nx( static_cast<int>(xDiv) );
        const int ny( static_cast<int>(yDiv) );

        const double xS( xDiv - nx );
        const double yS( yDiv - ny );

        const double c00( m_dataTable[ nx ][ ny ] );   // Lower left corner
        const double c01( m_dataTable[ nx ][ny+1] );   // Upper left corner
        const double c10( m_dataTable[nx+1][ ny ] );   // Lower right corner
        const double c11( m_dataTable[nx+1][ny+1] );   // Upper right corner

        return (c11 - c10 - c01 + c00)*xS*yS + (c10 - c00)*xS + (c01 - c00)*yS + c00;
    }

    double histogramI(const double x, const double y) const
    {
        if ( x < m_xStart or y < m_yStart or x > m_xStop or y > m_yStop )
            return 0.0;

        const double xDiv( (x - m_xMinI)/m_xStep );
        const double yDiv( (y - m_yMinI)/m_yStep );

        const int nx( (xDiv <= 0.0) ? 0 : (xDiv >= m_xRes-1) ? m_xRes-1 : static_cast<int>(xDiv) );
        const int ny( (yDiv <= 0.0) ? 0 : (yDiv >= m_yRes-1) ? m_yRes-1 : static_cast<int>(yDiv) );

        const double xS( (xDiv <= 0.0 or xDiv >= m_xRes-1) ? 0.0 : xDiv - nx );
        const double yS( (yDiv <= 0.0 or yDiv >= m_yRes-1) ? 0.0 : yDiv - ny );

        const double c00( m_dataTable[ nx ][ ny ] );   // Lower left corner
        const double c01( m_dataTable[ nx ][ny+1] );   // Upper left corner
        const double c10( m_dataTable[nx+1][ ny ] );   // Lower right corner
        const double c11( m_dataTable[nx+1][ny+1] );   // Upper right corner

        return (c11 - c10 - c01 + c00)*xS*yS + (c10 - c00)*xS + (c01 - c00)*yS + c00;
    }

    double operator() (const double x, const double y) const
    {
        return m_histogram ? histogramI(x, y) : linearI(x, y);
    }

    // Find the y value according to the 2D table using the cumulative distribution function
    double find_yFromCDF(const double x, double rnd) const
    {
        if ( x < m_xStart or x > m_xStop or rnd < 0.0 or rnd >= 1.0 )
            return 0;

        double ySumTable[m_yRes + 1];

        const double xDiv( (x - m_xMinI)/m_xStep );
        const int nx( (xDiv <= 0.0) ? 0 : (xDiv >= m_xRes - 1) ? m_xRes - 1 : static_cast<int>(xDiv) );
        const double xS( xDiv - nx );

        const double normalization( m_xDistributionTable[nx]*(1.0 - xS) + xS*m_xDistributionTable[nx + 1] );

        double sum(0.0), contrib(0.0), currentSum(0.0), previousSum(0.0);

        // Optimization by setting the values in ySumTable in one pass
        for (int kCnt(0); kCnt < m_yRes; ++kCnt)
        {
            contrib = ( m_dataTable[nx][kCnt]*(1.0 - xS) + xS*m_dataTable[nx + 1][kCnt] )*m_yStep/normalization;
            sum += contrib;
            currentSum = 0.5*sum; // Half now, half later
            ySumTable[kCnt] = previousSum + currentSum;
            previousSum = currentSum;
        }

        ySumTable[m_yRes] = sum; // There is no later

        // Bisection to find the index right below the rnd value
        int na( 0 );
        int nc( m_yRes );
        int nb( (na + nc)/2 );

        while (true)
        {
            if ( ySumTable[nb] - rnd < 0.0 )
            {
                if ( ySumTable[nb + 1] - rnd > 0.0 )
                    break;

                na = nb; // If the rnd number is larger, go to higher indices
            }
            else
                nc = nb;

            nb = (na + nc)/2;
        }

        const double y0( m_yStart + nb*m_yStep );
        const double u( ( rnd - ySumTable[nb] )/( ySumTable[nb + 1] - ySumTable[nb] ) ); // u is between 0 and 1, uniformly distributed

        if ( nb >= m_yRes - 1 )
            return y0 + m_yStep*0.5*u;

        const double z0( m_dataTable[nx][nb]*(1.0 - xS) + xS*m_dataTable[nx + 1][nb] );
        const double z1( m_dataTable[nx][nb + 1]*(1.0 - xS) + xS*m_dataTable[nx + 1][nb + 1] );
        const double dz( z1 - z0 );

        if ( std::abs(dz) < 1.0e-16 )
            return y0 + m_yStep*u;

        const double r( u*(z1 + z0)*dz ); // r is between 0 and (z0 + z1)*dz = z1^2 - z0^2, uniformly distributed

        // z0^2 + r is between z0^2 and z1^2, so it's always positive
        const double yS( r/dz/(sqrt( z0*z0 + r ) + z0) );

        const double y( y0 + m_yStep*yS );

        return y;
    }

    double* get_xDistribution() const
    {
        return m_xDistributionTable;
    }

    int get_xRes() const
    {
        return m_xRes;
    }

    double get_xStart() const
    {
        return m_xStart;
    }

    double get_xStep() const
    {
        return m_xStep;
    }

    double get_xStop() const
    {
        return m_xStop;
    }

    int get_yRes() const
    {
        return m_yRes;
    }

    double get_yStart() const
    {
        return m_yStart;
    }

    double get_yStep() const
    {
        return m_yStep;
    }

    double get_yStop() const
    {
        return m_yStop;
    }

    // Getter for normalization
    double get_normalization() const
    {
        return m_normalization;
    }

private:
    const bool m_histogram;
    bool m_valid; // Validity flag

    int m_xRes;
    double m_xStart;
    double m_xMinI; // Above this x, the values are interpolated
    double m_xStep;
    double m_xStop;

    int m_yRes;
    double m_yStart;
    double m_yMinI; // Above this y, the values are interpolated
    double m_yStep;
    double m_yStop;

    double m_normalization;

    double** m_dataTable;
    double* m_xDistributionTable;
};

#endif

