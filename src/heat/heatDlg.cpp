// heatDlg.cpp : implementation file
//

#include "stdafx.h"
#include "heat.h"
#include "heatDlg.h"
#include "afxdialogex.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// CHeatDlg dialog

using namespace plot;
using namespace util;
using namespace model;


CHeatDlg::CHeatDlg(CWnd* pParent /*=NULL*/)
    : CSimulationDialog(CHeatDlg::IDD, pParent)
    , m_cParams(make_default_parameters())
    , m_cPlotData(make_plot_data())
    , m_cIsotermsToDisplay(10)
    , m_cMaxTToDisplay(100)
    , m_cDisplayHeatMap(TRUE)
    , m_cDisplayHeatMapBool(true)
{
    m_hIcon = AfxGetApp()->LoadIcon(IDR_MAINFRAME);
}

void CHeatDlg::DoDataExchange(CDataExchange* pDX)
{
    CSimulationDialog::DoDataExchange(pDX);
    DDX_Control(pDX, IDC_PLOT, m_cPlot);
    DDX_Text(pDX, IDC_EDIT1, m_cParams.R);
    DDX_Text(pDX, IDC_EDIT2, m_cParams.L);
    DDX_Text(pDX, IDC_EDIT3, m_cParams.d);
    DDX_Text(pDX, IDC_EDIT4, m_cParams.R_h);
    DDX_Text(pDX, IDC_EDIT5, m_cParams.L_h);
    DDX_Text(pDX, IDC_EDIT12, m_cParams.z_h);
    DDX_Text(pDX, IDC_EDIT6, m_cParams.P_h);
    DDX_Text(pDX, IDC_EDIT7, m_cParams.d_c);
    DDX_Text(pDX, IDC_EDIT8, m_cParams.L_c);
    DDX_Text(pDX, IDC_EDIT9, m_cParams.z_c);
    DDX_Text(pDX, IDC_EDIT10, m_cParams.dt);
    DDX_Text(pDX, IDC_EDIT11, m_cParams.lambda_m);
    DDX_Text(pDX, IDC_EDIT13, m_cIsotermsToDisplay);
    DDX_Text(pDX, IDC_EDIT14, m_cMaxTToDisplay);
    DDX_Check(pDX, IDC_CHECK1, m_cDisplayHeatMap);
}

BEGIN_MESSAGE_MAP(CHeatDlg, CSimulationDialog)
    ON_WM_PAINT()
    ON_WM_QUERYDRAGICON()
    ON_BN_CLICKED(IDC_BUTTON1, &CHeatDlg::OnBnClickedButton1)
    ON_BN_CLICKED(IDC_BUTTON2, &CHeatDlg::OnBnClickedButton2)
    ON_BN_CLICKED(IDC_CHECK1, &CHeatDlg::OnBnClickedCheck1)
END_MESSAGE_MAP()


// CHeatDlg message handlers

BOOL CHeatDlg::OnInitDialog()
{
    CSimulationDialog::OnInitDialog();

    // Set the icon for this dialog.  The framework does this automatically
    //  when the application's main window is not a dialog
    SetIcon(m_hIcon, TRUE);            // Set big icon
    SetIcon(m_hIcon, FALSE);        // Set small icon

    // TODO: Add extra initialization here

    adjust(m_cParams, m_cPlotData);
    m_cPlot.plot_layer.with(make_root_drawable(m_cPlotData, {
        custom_drawable::create(make_system_painter(m_cParams))
    }));

    m_cPlot.background = palette::brush();

    m_cPlot.RedrawWindow();

    return TRUE;  // return TRUE  unless you set the focus to a control
}

// If you add a minimize button to your dialog, you will need the code below
//  to draw the icon.  For MFC applications using the document/view model,
//  this is automatically done for you by the framework.

void CHeatDlg::OnPaint()
{
    if (IsIconic())
    {
        CPaintDC dc(this); // device context for painting

        SendMessage(WM_ICONERASEBKGND, reinterpret_cast<WPARAM>(dc.GetSafeHdc()), 0);

        // Center icon in client rectangle
        int cxIcon = GetSystemMetrics(SM_CXICON);
        int cyIcon = GetSystemMetrics(SM_CYICON);
        CRect rect;
        GetClientRect(&rect);
        int x = (rect.Width() - cxIcon + 1) / 2;
        int y = (rect.Height() - cyIcon + 1) / 2;

        // Draw the icon
        dc.DrawIcon(x, y, m_hIcon);
    }
    else
    {
        CSimulationDialog::OnPaint();
    }
}

// The system calls this function to obtain the cursor to display while the user drags
//  the minimized window.
HCURSOR CHeatDlg::OnQueryDragIcon()
{
    return static_cast<HCURSOR>(m_hIcon);
}


void CHeatDlg::OnBnClickedButton1()
{
    UpdateData(TRUE);

    StartSimulationThread();
}


void CHeatDlg::OnBnClickedButton2()
{
    StopSimulationThread();
}

void CHeatDlg::OnBnClickedCheck1()
{
    UpdateData(TRUE);
    this->m_cDisplayHeatMapBool = (this->m_cDisplayHeatMap != 0);
}
