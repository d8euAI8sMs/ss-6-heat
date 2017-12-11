// heatDlg.h : header file
//

#pragma once

#include <util/common/gui/SimulationDialog.h>
#include <util/common/plot/PlotStatic.h>

#include "model.h"

// CHeatDlg dialog
class CHeatDlg : public CSimulationDialog
{
// Construction
public:
    CHeatDlg(CWnd* pParent = NULL);    // standard constructor

// Dialog Data
    enum { IDD = IDD_HEAT_DIALOG };

    protected:
    virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support


// Implementation
protected:
    HICON m_hIcon;

    // Generated message map functions
    virtual BOOL OnInitDialog();
    afx_msg void OnPaint();
    afx_msg HCURSOR OnQueryDragIcon();
    DECLARE_MESSAGE_MAP()
public:
    virtual void OnSimulation();
    afx_msg void OnBnClickedButton1();
    afx_msg void OnBnClickedButton2();
    PlotStatic m_cPlot;
    model::parameters m_cParams;
    model::plot_data m_cPlotData;
    model::chasing_data m_cChasingData;
    std::vector < std::vector < double > > T;
    int m_cIsotermsToDisplay;
    double m_cMaxTToDisplay;
    double m_cdT;
    BOOL m_cDisplayHeatMap;
    bool m_cDisplayHeatMapBool;
    afx_msg void OnBnClickedCheck1();
};
