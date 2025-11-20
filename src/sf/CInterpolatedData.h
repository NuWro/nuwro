#ifndef C_INTERPOLATED_DATA_H
#define C_INTERPOLATED_DATA_H

#include <fstream>
#include <iostream>
#include <cstdlib>

///   Class handling 1D data tables

///   The input format is the following:
///
///   first the header:
///   xRes xStep
///
///   then xRes rows of data:
///   the x value and the corresponding function value

class CInterpolatedData
{
public:
   CInterpolatedData(const std::string inputData, const double i_abscissaUnit, const double i_ordinateUnit, const bool i_histogram)
   :m_abscissaUnit(i_abscissaUnit),
    m_ordinateUnit(i_ordinateUnit),
    m_histogram(i_histogram)
   {
        std::ifstream dataFile( inputData.c_str() );

        if (dataFile.fail())
        {
            std::cerr<<"Indispensable file '"<<inputData<<"' not found"<<std::endl;
          //  system("pause");
            return;
        }

        double dummy(0.0);
        dataFile>>m_xRes>>m_xStep;
        m_xStep *= m_abscissaUnit;
        m_dataTable = new double[m_xRes];

        for (int lCnt(0); lCnt < m_xRes; ++lCnt)
        {
            dataFile>>dummy;
            if (lCnt == 0)
            {
                m_xMinI = dummy*m_abscissaUnit;
                m_xStart = m_histogram ? m_xMinI - 0.5*m_xStep : m_xMinI;
            }
            else
                if (lCnt == (m_xRes - 1))
                {
                    m_xMaxI = dummy*m_abscissaUnit;
                    m_xStop = m_histogram ? m_xMaxI + 0.5*m_xStep : m_xMaxI;
                }

            dataFile>>dummy;
            m_dataTable[lCnt] = dummy*m_ordinateUnit;
        }
    //  std::cout<<"The file '"<<inputData<<"' is loaded"<<std::endl;
        if ( m_histogram )
            m_eval = &CInterpolatedData::histogramI;
        else
            m_eval = &CInterpolatedData::linearI;
    };

    CInterpolatedData(const double* i_inputData, const int i_xRes, const double i_abscissaUnit, const double i_ordinateUnit, const double i_xStart, const double i_xStep, const double i_xStop, const bool i_histogram)
   :m_abscissaUnit(i_abscissaUnit),
    m_ordinateUnit(i_ordinateUnit),
    m_histogram(i_histogram),
    m_xRes(i_xRes),
    m_xMinI(i_xStart*i_abscissaUnit),
    m_xStart(m_histogram ? m_xMinI - 0.5*i_xStep*i_abscissaUnit : m_xMinI),
    m_xStep(i_xStep*i_abscissaUnit),
    m_xMaxI(i_xStop*i_abscissaUnit),
    m_xStop(m_histogram ? m_xMaxI + 0.5*i_xStep*i_abscissaUnit : m_xMaxI)
   {
        m_dataTable = new double[m_xRes];

        for ( int lCnt(0); lCnt < m_xRes; ++lCnt )
            m_dataTable[lCnt] = i_inputData[lCnt]*m_ordinateUnit;

        if ( m_histogram )
            m_eval = &CInterpolatedData::histogramI;
        else
            m_eval = &CInterpolatedData::linearI;
    };

    ~CInterpolatedData()
    {
        delete []m_dataTable;
    };

    typedef double (CInterpolatedData::*DEvalPtr)(const double) const;

    double withoutI(const double x) const
    {
        if ( x < m_xStart or x > m_xStop )
            return 0.0;

        const double xDiv( (x - m_xStart)/m_xStep );
        const int nx( static_cast<int>(xDiv) );
        const double c0( m_dataTable[nx] );
        return c0;
    }

    double linearI(const double x) const
    {
        if ( x >= m_xMaxI )
            return m_dataTable[m_xRes - 1];
        if ( x <= m_xMinI )
            return m_dataTable[0];

        const double xDiv( (x - m_xMinI)/m_xStep );
        const int nx( static_cast<int>(xDiv) );
        const double xS( xDiv - nx );
        const double c0( m_dataTable[nx] );   ///lower corner
        const double c1( m_dataTable[nx+1] ); ///upper corner
        return (c1 - c0)*xS + c0;
    }

    double histogramI(const double x) const
    {
        if ( x >= m_xMaxI )
        {
            if ( x > m_xStop )
                return 0.0;
            return m_dataTable[m_xRes - 1];
        }

        if ( x <= m_xMinI )
        {
            if ( x < m_xStart )
                return 0.0;
            return m_dataTable[0];
        }

        const double xDiv( (x - m_xMinI)/m_xStep );
        const int nx( static_cast<int>(xDiv) );
        const double xS( xDiv - nx );
        const double c0( m_dataTable[nx] );   ///lower corner
        const double c1( m_dataTable[nx+1] ); ///upper corner
        return (c1 - c0)*xS + c0;
    }

    DEvalPtr m_eval;
    double eval(const double x) const
    {
       return (this->*m_eval) (x);
    }

    double operator() (const double x) const
    {
        return eval(x);
    }

    int get_xRes() const
    {
        return m_xRes;
    }

    double get_xMinI() const
    {
        return m_xMinI;
    }

    double get_xStart() const
    {
        return m_xStart;
    }

    double get_xStep() const
    {
        return m_xStep;
    }

    double get_xMaxI() const
    {
        return m_xMaxI;
    }

    double get_xStop() const
    {
        return m_xStop;
    }

    double get_yStart() const
    {
        return m_dataTable[0];
    }

    double get_yStop() const
    {
        return m_dataTable[m_xRes-1];
    }

private:
    const double m_abscissaUnit;
    const double m_ordinateUnit;
    const bool m_histogram;

    mutable int m_xRes;
    mutable double m_xMinI;///above this x, the values are interpolated
    mutable double m_xStart;
    mutable double m_xStep;
    mutable double m_xMaxI;///below this x, the values are interpolated
    mutable double m_xStop;
    mutable double* m_dataTable;
};
#endif
