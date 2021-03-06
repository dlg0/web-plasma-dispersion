#include <Wt/WApplication>
#include <Wt/WBreak>
#include <Wt/WContainerWidget>
#include <Wt/WLineEdit>
#include <Wt/WPushButton>
#include <Wt/WText>
#include <Wt/WDoubleValidator>
#include <Wt/WValidator>
#include <Wt/WColor>
#include <Wt/WComboBox>
#include <Wt/WGroupBox>
#include <Wt/WButtonGroup>
#include <Wt/WRadioButton>
#include <sstream>
#include <Wt/WBorderLayout>
#include <Wt/WImage>
#include <Wt/Chart/WCartesianChart>
#include <Wt/WStandardItemModel>
#include <boost/lexical_cast.hpp>
#include <Wt/Chart/WAxis>
#include <stixVars.hpp>
#include <constants.hpp>
#include "dielectric.hpp"
#include <armadillo>
#include "rotation.hpp"
#include <o2scl/poly.h>
//#include "zcircle.hpp"
#define GNUPLOT_ENABLE_BLITZ
#include "gnuplot-iostream.h"
//#include <boost/iostreams/device/file_descriptor.hpp>
#include <blitz/array.h>
#include <map>
#include <cstdio>
#include <assert.h>

class IonSpec
{
		public:
		IonSpec(int, Wt::WGroupBox*);
		~IonSpec(void);
		void setVisible(void);
		void setHidden(void);

		private:

		public:
		Wt::WLineEdit *SpecAMU_LineEdit, *SpecZ_LineEdit, *SpecT_eV_LineEdit, 
				*SpecDensity_m3_LineEdit, *SpecHarmonicN_LineEdit;
		Wt::WBreak *SpecBreak;
		Wt::WText *SpecText;
		int SpecNo;
};

IonSpec::IonSpec(int SpecNoIn, Wt::WGroupBox *IonSpecies_Container )
{
	SpecNo = SpecNoIn;
	std::stringstream TmpStr;
	TmpStr << "Ion Species " << SpecNoIn << " (AMU, Z, T [eV], Density [m^-3], MaxHarmN)";
	Wt::WString TmpStrWt(TmpStr.str());
	SpecText = new Wt::WText(TmpStrWt,IonSpecies_Container);
	SpecText->setHidden(1);
	SpecAMU_LineEdit = new Wt::WLineEdit("1",IonSpecies_Container);
	SpecAMU_LineEdit->setHidden(1);
	SpecZ_LineEdit = new Wt::WLineEdit("1",IonSpecies_Container);
	SpecZ_LineEdit->setHidden(1);
	SpecT_eV_LineEdit = new Wt::WLineEdit("0.01",IonSpecies_Container);
	SpecT_eV_LineEdit->setHidden(1);
	SpecDensity_m3_LineEdit = new Wt::WLineEdit("1e19",IonSpecies_Container);
	SpecDensity_m3_LineEdit->setHidden(1);
	SpecHarmonicN_LineEdit = new Wt::WLineEdit("3",IonSpecies_Container);
	SpecHarmonicN_LineEdit->setHidden(1);
	SpecBreak = new Wt::WBreak(IonSpecies_Container);
SpecBreak->setHidden(1);
}
IonSpec::~IonSpec(void)
{
}
void IonSpec::setVisible(void)
{
	SpecText->setHidden(0);
	SpecAMU_LineEdit->setHidden(0);
	SpecZ_LineEdit->setHidden(0);
	SpecT_eV_LineEdit->setHidden(0);
	SpecDensity_m3_LineEdit->setHidden(0);
	SpecHarmonicN_LineEdit->setHidden(0);
	SpecBreak->setHidden(0);
}
void IonSpec::setHidden(void)
{
	SpecText->setHidden(1);
	SpecAMU_LineEdit->setHidden(1);
	SpecZ_LineEdit->setHidden(1);
	SpecT_eV_LineEdit->setHidden(1);
	SpecDensity_m3_LineEdit->setHidden(1);
	SpecHarmonicN_LineEdit->setHidden(1);
	SpecBreak->setHidden(1);
}

class DispersionApp : public Wt::WApplication
{
		public:
		    DispersionApp(const Wt::WEnvironment& env);

		private:
			Wt::WPushButton *Update_PushButton;
			Wt::WLineEdit *nameEdit_;
			Wt::WText *greeting_;
			Wt::WLineEdit *BField_LineEdit, *Freq_Hz_LineEdit;
			Wt::WDoubleValidator *BFieldValidator, *FreqValidator;
			Wt::WText *BFieldValidText, *FreqValidText;

			Wt::WContainerWidget *w;
			Wt::WBorderLayout *layout;

			std::vector<IonSpec> IonSpecies;

			Wt::WGroupBox *West_Container, *Center_Container, *East_Container, *IonSpecies_Container, 
					*PointParameterType_Container,*MagneticParameterType_Container;
			Wt::WButtonGroup *ParameterType_ButtonGroup, *yRange_ButtonGroup;
			Wt::WComboBox *nIonSpecies_ComboBox;

			// Scan
			int nX;
			Wt::WText *nX_Text;
			Wt::WLineEdit *nX_LineEdit;

			// Magnetic field scan
			float r0, b0, rMin, rMax, kPar, kyRe,kyIm,kzRe,kzIm, yRangeMin, yRangeMax, bxFrac, byFrac;
			std::complex<double> ky, kz;
			Wt::WText *r0_Text, *b0_Text, *rMin_Text, *rMax_Text, *kPar_Text,
					*ky_Text, *kz_Text, *bxFrac_Text, *byFrac_Text;
			Wt::WLineEdit *r0_LineEdit, *b0_LineEdit, *rMin_LineEdit, *rMax_LineEdit, 
					*kPar_LineEdit, *kzRe_LineEdit, *kyRe_LineEdit, *bxFrac_LineEdit, *byFrac_LineEdit,
					*kzIm_LineEdit, *kyIm_LineEdit;

			// electron parameters
			Wt::WLineEdit *electron_AMU_LineEdit, *electron_Z_LineEdit, 
					*electron_T_eV_LineEdit, *electron_Density_m3_LineEdit,
					*electron_HarmonicN_LineEdit;

			Wt::WLineEdit *nuOmg_LineEdit;

			// Plot Y-Range controllers
			Wt::WText *yRangeMin_Text, *yRangeMax_Text;
			Wt::WLineEdit *yRangeMin_LineEdit, *yRangeMax_LineEdit;

			// k-space scan parameters
			Wt::WText *kxReMin_Text, *kxReMax_Text, *kxImMin_Text, *kxImMax_Text, 
					*kxReN_Text, *kxImN_Text;
			Wt::WLineEdit *kxReMin_LineEdit, *kxReMax_LineEdit, *kxImMin_LineEdit, *kxImMax_LineEdit,
					*kxReN_LineEdit, *kxImN_LineEdit;

			std::vector <float> x;//,b;
			arma::colvec bx,by,bz,bMag;

			Wt::Chart::WCartesianChart *b_CartesianChart, *a_CartesianChart;
			Wt::WStandardItemModel *b_model, *a_model;
			Wt::Chart::WDataSeries *b_series; 
			Wt::Chart::WDataSeries *a_series1,*a_series2,*a_series3,*a_series4;
			Wt::Chart::WDataSeries *a_series5,*a_series6,*a_series7,*a_series8;

			Wt::Chart::WAxis *aYAxis, *aXAxis, *bYAxis, *bXAxis;

			int nSpecies;

			enum yRangeType { Auto = 0, Manual = 1 };

			void greet();
			void editedLineEdit(Wt::WLineEdit *LineEditIn, Wt::WText *ValidText);
			void UpdateIonSpecies();
			void RevealParameterSettings();
			void UpdateCalculation();
			void UpdateBField();
			void SetYRange();
			void AddChart(std::string title, std::string xTitle, std::string yTitle, 
							Wt::Chart::WCartesianChart *_chart, Wt::WStandardItemModel *_model, 
							Wt::Chart::WAxis *_xAxis, Wt::Chart::WAxis *_yAxis);
			void SetYRangeManual();


};

DispersionApp::DispersionApp(const Wt::WEnvironment& env) : Wt::WApplication(env)
{
		    setTitle("Plasma Dispersion Calculator");

			w = new Wt::WContainerWidget(root());
			layout = new Wt::WBorderLayout();
			w->setLayout(layout,Wt::AlignTop | Wt::AlignJustify);

			//layout->addWidget(new Wt::WImage("naruto.jpg"),Wt::WBorderLayout::Center);

			enum CalculationType { Hot = 1, Cold = 2 };

			West_Container = new Wt::WGroupBox("West");
			West_Container->setId("West");
			West_Container->resize(300,400);
			layout->addWidget(West_Container,Wt::WBorderLayout::West);

			East_Container = new Wt::WGroupBox("East");
			East_Container->setId("East");
			East_Container->resize(250,400);
			layout->addWidget(East_Container,Wt::WBorderLayout::East);

			Center_Container = new Wt::WGroupBox("Center");
			Center_Container->setId("Center");
			layout->addWidget(Center_Container,Wt::WBorderLayout::Center);

			Wt::WButtonGroup *HotOrCold_ButtonGroup = new Wt::WButtonGroup(West_Container);
			Wt::WRadioButton *HotOrCold_Button;

			HotOrCold_Button = new Wt::WRadioButton("Hot", West_Container);
			new Wt::WBreak(West_Container);
			HotOrCold_ButtonGroup->addButton(HotOrCold_Button,Hot);

			HotOrCold_Button = new Wt::WRadioButton("Cold", West_Container);
			new Wt::WBreak(West_Container);
			HotOrCold_ButtonGroup->addButton(HotOrCold_Button,Cold);

			HotOrCold_ButtonGroup->setCheckedButton(HotOrCold_ButtonGroup->button(Cold));

			nX_Text = new Wt::WText("nX",West_Container);
			nX_LineEdit = new Wt::WLineEdit("100",West_Container);
			nX = boost::lexical_cast<int>(nX_LineEdit->text().narrow());

			Wt::WRadioButton *yRangeAuto_Button, *yRangeManu_Button;
			yRange_ButtonGroup = new Wt::WButtonGroup(East_Container);
		   	yRangeAuto_Button = new Wt::WRadioButton("Auto y-range",East_Container);
			yRange_ButtonGroup->addButton(yRangeAuto_Button);
			new Wt::WBreak(East_Container);
			yRangeManu_Button = new Wt::WRadioButton("Manual y-range",East_Container);
			yRange_ButtonGroup->addButton(yRangeManu_Button);
			new Wt::WBreak(East_Container);
			yRange_ButtonGroup->setCheckedButton(yRange_ButtonGroup->button(Auto));
			yRange_ButtonGroup->checkedChanged().connect(
							boost::bind(&DispersionApp::SetYRange,this));


			yRangeMin_Text = new Wt::WText("yRangeMin",East_Container);
			yRangeMin_LineEdit = new Wt::WLineEdit("-5000",East_Container);
			yRangeMin_LineEdit->changed().connect(
							boost::bind(&DispersionApp::SetYRange,this));
			yRangeMin = boost::lexical_cast<double>(yRangeMin_LineEdit->text().narrow());
			new Wt::WBreak(East_Container);

			yRangeMax_Text = new Wt::WText("yRangeMax",East_Container);
			yRangeMax_LineEdit = new Wt::WLineEdit("5000",East_Container);
			yRangeMax_LineEdit->changed().connect(
							boost::bind(&DispersionApp::SetYRange,this));
			yRangeMax = boost::lexical_cast<double>(yRangeMax_LineEdit->text().narrow());
			new Wt::WBreak(East_Container);

			enum ParameterType { Point = 1, MagneticField = 2, Numerical = 3 };

			//ParameterType_Container = new Wt::WGroupBox("Parameter Scan Type");
			//layout->addWidget(ParameterType_Container,Wt::WBorderLayout::East);
			ParameterType_ButtonGroup = new Wt::WButtonGroup(East_Container);

			Wt::WRadioButton *ParameterType_Button;

			ParameterType_Button = new Wt::WRadioButton("Point", East_Container);
			new Wt::WBreak(East_Container);
			ParameterType_ButtonGroup->addButton(ParameterType_Button,Point);

			ParameterType_Button = new Wt::WRadioButton("MagneticField", East_Container);
			new Wt::WBreak(East_Container);
			ParameterType_ButtonGroup->addButton(ParameterType_Button,MagneticField);

			ParameterType_Button = new Wt::WRadioButton("Numerical", East_Container);
			new Wt::WBreak(East_Container);
			ParameterType_ButtonGroup->addButton(ParameterType_Button,Numerical);

			kxReN_Text = new Wt::WText("n Re(kx):",East_Container);
			kxReN_LineEdit = new Wt::WLineEdit("200",East_Container);
			new Wt::WBreak(East_Container);

			kxImN_Text = new Wt::WText("n Im(kx):",East_Container);
			kxImN_LineEdit = new Wt::WLineEdit("100",East_Container);
			new Wt::WBreak(East_Container);

			kxReMin_Text = new Wt::WText("Re(kx) min:",East_Container);
			kxReMin_LineEdit = new Wt::WLineEdit("-40",East_Container);
			new Wt::WBreak(East_Container);

			kxReMax_Text = new Wt::WText("Re(kx) Max:",East_Container);
			kxReMax_LineEdit = new Wt::WLineEdit("40",East_Container);
			new Wt::WBreak(East_Container);

			kxImMin_Text = new Wt::WText("Im(kx) min:",East_Container);
			kxImMin_LineEdit = new Wt::WLineEdit("-20",East_Container);
			new Wt::WBreak(East_Container);

			kxImMax_Text = new Wt::WText("Im(kx) Max:",East_Container);
			kxImMax_LineEdit = new Wt::WLineEdit("20",East_Container);
			new Wt::WBreak(East_Container);


			ParameterType_ButtonGroup->setCheckedButton(ParameterType_ButtonGroup->button(MagneticField));
			ParameterType_ButtonGroup->checkedChanged().connect(
							boost::bind(&DispersionApp::RevealParameterSettings,this));

			PointParameterType_Container = new Wt::WGroupBox("Point Parameter Settings",root());

			new Wt::WBreak(West_Container);

			Wt::WText *BFieldText = new Wt::WText(PointParameterType_Container);
			BFieldText->setText("BField [T]");
			BFieldValidator = new Wt::WDoubleValidator(1e-5,5e5,PointParameterType_Container);
			BFieldValidator->setMandatory(1);
			BField_LineEdit = new Wt::WLineEdit("5",PointParameterType_Container);
			BField_LineEdit->setValidator(BFieldValidator);
			BFieldValidText = new Wt::WText(PointParameterType_Container);
			BField_LineEdit->changed().connect(
					boost::bind(&DispersionApp::editedLineEdit,this,BField_LineEdit,BFieldValidText));
			//new Wt::WBreak(West_Container);

			Wt::WText *FreqText = new Wt::WText(West_Container);
			FreqText->setText("Freq [Hz]");
			FreqValidator = new Wt::WDoubleValidator(1e1,5e12,West_Container);
			FreqValidator->setMandatory(1);
			Freq_Hz_LineEdit = new Wt::WLineEdit("30e6",West_Container);
			Freq_Hz_LineEdit->setValidator(FreqValidator);
			FreqValidText = new Wt::WText(West_Container);
			Freq_Hz_LineEdit->changed().connect(
					boost::bind(&DispersionApp::editedLineEdit,this,Freq_Hz_LineEdit,FreqValidText));
			new Wt::WBreak(West_Container);
			Wt::WText *nuOmgText = new Wt::WText("nu/Omega",West_Container);
			nuOmg_LineEdit = new Wt::WLineEdit("0.02",West_Container);
			new Wt::WBreak(West_Container);

			PointParameterType_Container->setHidden(1);

			MagneticParameterType_Container = new Wt::WGroupBox("Magnetic Field Scan Settings",root());
			MagneticParameterType_Container->setHidden(0);

			//kPar_Text = new Wt::WText("kPar",West_Container);
			//kPar_LineEdit = new Wt::WLineEdit("0.1",West_Container);
			//new Wt::WBreak(West_Container);

			ky_Text = new Wt::WText("ky[Re,Im]",West_Container);
			kyRe_LineEdit = new Wt::WLineEdit("0",West_Container);
			kyIm_LineEdit = new Wt::WLineEdit("0",West_Container);
			new Wt::WBreak(West_Container);

			kz_Text = new Wt::WText("kz[Re,Im]",West_Container);
			kzRe_LineEdit = new Wt::WLineEdit("10",West_Container);
			kzIm_LineEdit = new Wt::WLineEdit("0",West_Container);
			new Wt::WBreak(West_Container);

			r0_Text = new Wt::WText("r0",West_Container);
			r0_LineEdit = new Wt::WLineEdit("2",West_Container);
			new Wt::WBreak(West_Container);

			b0_Text = new Wt::WText("b0",West_Container);
			b0_LineEdit = new Wt::WLineEdit("5",West_Container);
			new Wt::WBreak(West_Container);

			bxFrac_Text = new Wt::WText("bxFrac",West_Container);
			bxFrac_LineEdit = new Wt::WLineEdit("0.0",West_Container);
			new Wt::WBreak(West_Container);

			byFrac_Text = new Wt::WText("byFrac",West_Container);
			byFrac_LineEdit = new Wt::WLineEdit("0.0",West_Container);
			new Wt::WBreak(West_Container);

			rMin_Text = new Wt::WText("rMin",West_Container);
			rMin_LineEdit = new Wt::WLineEdit("0.1",West_Container);
			new Wt::WBreak(West_Container);

			rMax_Text = new Wt::WText("rMax",West_Container);
			rMax_LineEdit = new Wt::WLineEdit("4",West_Container);

			Update_PushButton = new Wt::WPushButton("Update",West_Container);
			Update_PushButton->clicked().connect(this, &DispersionApp::UpdateCalculation);

			IonSpecies_Container = new Wt::WGroupBox("Ion Species",root());

			Wt::WText *electronText = new Wt::WText("Electron parameters (AMU, Z, T [eV], Density [m^-3], MaxHarmN)",IonSpecies_Container);
			electron_AMU_LineEdit = new Wt::WLineEdit("0.000544617",IonSpecies_Container);
			electron_Z_LineEdit = new Wt::WLineEdit("-1",IonSpecies_Container);
			electron_T_eV_LineEdit = new Wt::WLineEdit("0.01",IonSpecies_Container);
			electron_Density_m3_LineEdit = new Wt::WLineEdit("auto",IonSpecies_Container);
			electron_HarmonicN_LineEdit = new Wt::WLineEdit("3",IonSpecies_Container);

			new Wt::WBreak(IonSpecies_Container);

			Wt::WText *nIonSpecText = new Wt::WText("N Ion Species",IonSpecies_Container);
			nIonSpecies_ComboBox = new Wt::WComboBox(IonSpecies_Container);
			nIonSpecies_ComboBox->addItem("1");
			nIonSpecies_ComboBox->addItem("2");
			nIonSpecies_ComboBox->addItem("3");
			nIonSpecies_ComboBox->addItem("4");
			nIonSpecies_ComboBox->addItem("5");
			nIonSpecies_ComboBox->addItem("6");
			nIonSpecies_ComboBox->setCurrentIndex(0);
			nIonSpecies_ComboBox->changed().connect(
				boost::bind(&DispersionApp::UpdateIonSpecies,this));

			new Wt::WBreak(IonSpecies_Container);


			a_CartesianChart = new Wt::Chart::WCartesianChart(Center_Container); 
			b_CartesianChart = new Wt::Chart::WCartesianChart(Center_Container); 

			a_model = new Wt::WStandardItemModel(nX,9);
			b_model = new Wt::WStandardItemModel(nX,2);

			a_series1 = new Wt::Chart::WDataSeries(1,Wt::Chart::PointSeries);
			a_series2 = new Wt::Chart::WDataSeries(2,Wt::Chart::PointSeries);
			a_series3 = new Wt::Chart::WDataSeries(3,Wt::Chart::PointSeries);
			a_series4 = new Wt::Chart::WDataSeries(4,Wt::Chart::PointSeries);
			a_series5 = new Wt::Chart::WDataSeries(5,Wt::Chart::PointSeries);
			a_series6 = new Wt::Chart::WDataSeries(6,Wt::Chart::PointSeries);
			a_series7 = new Wt::Chart::WDataSeries(7,Wt::Chart::PointSeries);
			a_series8 = new Wt::Chart::WDataSeries(8,Wt::Chart::PointSeries);

			b_series = new Wt::Chart::WDataSeries(1,Wt::Chart::LineSeries);

			aYAxis = &(a_CartesianChart->axis(Wt::Chart::YAxis)); 
			aXAxis = &(a_CartesianChart->axis(Wt::Chart::XAxis)); 

			bYAxis = &(b_CartesianChart->axis(Wt::Chart::YAxis)); 
			bXAxis = &(b_CartesianChart->axis(Wt::Chart::XAxis)); 

			AddChart(std::string("a"),std::string("x [m]"),std::string("b [T]"),
						   a_CartesianChart, a_model, aXAxis, aYAxis);
			AddChart(std::string("b"),std::string("x [m]"),std::string("b [T]"),
						   b_CartesianChart, b_model, bXAxis, bYAxis);

			a_CartesianChart->addSeries(*a_series1);
			a_CartesianChart->addSeries(*a_series2);
			a_CartesianChart->addSeries(*a_series3);
			a_CartesianChart->addSeries(*a_series4);
			a_CartesianChart->addSeries(*a_series5);
			a_CartesianChart->addSeries(*a_series6);
			a_CartesianChart->addSeries(*a_series7);
			a_CartesianChart->addSeries(*a_series8);

			b_CartesianChart->addSeries(*b_series);

			nSpecies = 0;
			UpdateIonSpecies();
			UpdateBField();
			//UpdateCalculation();
}

void DispersionApp::AddChart (std::string title, std::string xTitle, std::string yTitle, 
	 Wt::Chart::WCartesianChart *_Chart, Wt::WStandardItemModel *_model, 
		Wt::Chart::WAxis *_xAxis, Wt::Chart::WAxis *_yAxis )
{

	//int nX = x.size();
	//_model->setHeaderData(0, Wt::WString("x (m)"));
	//_model->setHeaderData(1, Wt::WString("b (T)"));

	for (unsigned i = 0; i < nX; ++i) {
	  _model->setData(i, 0, i);
	  _model->setData(i, 1, i);
	}

	_Chart->setModel(_model);
	_Chart->setType(Wt::Chart::ScatterPlot);
	//_Chart->setLegendEnabled(true);
	_Chart->axis(Wt::Chart::XAxis).setLocation(Wt::Chart::ZeroValue);
	_Chart->resize(450,300);
	_Chart->setPlotAreaPadding(100,Wt::Left);
	_Chart->setPlotAreaPadding(40,Wt::Bottom);
	_Chart->setLegendLocation(Wt::Chart::LegendOutside,Wt::Top,Wt::AlignRight);
	_Chart->setTitle(Wt::WString(title));

	//Wt::Chart::WAxis &yAxis = _Chart->axis(Wt::Chart::YAxis); 
	//Wt::Chart::WAxis &xAxis = _Chart->axis(Wt::Chart::XAxis); 

	//_yAxis = _Chart->axis(Wt::Chart::YAxis); 
	//_xAxis = _Chart->axis(Wt::Chart::XAxis); 
	
	_yAxis->setTitle(yTitle);
	_xAxis->setTitle(xTitle);

	_Chart->setXSeriesColumn(0);
	//_Chart->addSeries(s);
}

void DispersionApp::RevealParameterSettings()
{
		if(ParameterType_ButtonGroup->checkedId()==1)
		{
			PointParameterType_Container->setHidden(0);
		}

		if(ParameterType_ButtonGroup->checkedId()==2)
		{
			PointParameterType_Container->setHidden(1);
		}

		if(ParameterType_ButtonGroup->checkedId()==3)
		{
			PointParameterType_Container->setHidden(1);
		}
}

void DispersionApp::UpdateIonSpecies()
{

		int nSpecTmp = nIonSpecies_ComboBox->currentIndex();
		std::cout << nSpecTmp << std::endl;

		for(int s=0; s<nSpecies; s++)
		{
				IonSpecies[s].setHidden();
		}

		IonSpecies.clear();
		nSpecies = 0;

		for(int s=0; s<=nSpecTmp; s++)
		{
				IonSpec mySpecParams(s+1,IonSpecies_Container);	
				IonSpecies.push_back(mySpecParams);
				IonSpecies[s].setVisible();
				nSpecies++;
		}

}

void DispersionApp::editedLineEdit(Wt::WLineEdit *LineEditIn, Wt::WText *ValidText)
{
		if(LineEditIn->validate() != Wt::WValidator::Valid)
		{
				LineEditIn->decorationStyle().setBackgroundColor(Wt::WColor("#F79A9A"));
				ValidText->setText("Not Valid Input");
		}
		else
		{
			LineEditIn->decorationStyle().setBackgroundColor(Wt::WColor("#66FF66"));
			ValidText->setText("Good to go");
		}
}

void DispersionApp::greet()
{
		    greeting_->setText("Hello there, " + nameEdit_->text());
}

void DispersionApp::UpdateBField ()
{
			nX = boost::lexical_cast<int>(nX_LineEdit->text().narrow());
			rMin = boost::lexical_cast<float>(rMin_LineEdit->text().narrow());
			rMax = boost::lexical_cast<float>(rMax_LineEdit->text().narrow());
			b0 = boost::lexical_cast<float>(b0_LineEdit->text().narrow());
			bxFrac = boost::lexical_cast<float>(bxFrac_LineEdit->text().narrow());
			byFrac = boost::lexical_cast<float>(byFrac_LineEdit->text().narrow());

			try {
				bxFrac = boost::lexical_cast<float>(bxFrac_LineEdit->text().narrow());
				std::cout << "Success: bxFrac_LineEdit read as " << bxFrac << std::endl;
			} catch(boost::bad_lexical_cast &) {
				std::cout << "Error: bxFrac_LineEdit could not be read" << std::endl;
			}

			r0 = boost::lexical_cast<float>(r0_LineEdit->text().narrow());

			x.resize(nX);
			//b.resize(nX);
			bx.set_size(nX);
			by.set_size(nX);
			bz.set_size(nX);
			bMag.set_size(nX);

			for(int i=0;i<nX;i++)
			{
					x[i] = i*(rMax-rMin)/nX+rMin;
					//b[i] = b0*r0/x[i];
					bz(i) = b0*r0/x[i];
					bx(i) = bxFrac*bz(i);
					by(i) = byFrac*bz(i);
					bMag(i) = sqrt(pow(bx(i),2)+pow(by(i),2)+pow(bz(i),2));

					assert(bMag(i)==bMag(i));
			}
}

void DispersionApp::SetYRangeManual()
{
		yRange_ButtonGroup->setCheckedButton(yRange_ButtonGroup->button(Manual));
}

void DispersionApp::SetYRange()
{
	
		int _auto = yRange_ButtonGroup->checkedId();
		std::cout << "yRange_ButtonGroup Id: " << _auto << std::endl;

		if(_auto!=0)
		{
			yRangeMin = boost::lexical_cast<double>(yRangeMin_LineEdit->text().narrow());
			yRangeMax = boost::lexical_cast<double>(yRangeMax_LineEdit->text().narrow());

			aYAxis->setMaximum(yRangeMax);
			aYAxis->setMinimum(yRangeMin);
		}
		else
		{
			aYAxis->setAutoLimits(Wt::Chart::MinimumValue | Wt::Chart::MaximumValue);
		}
}

std::vector<std::complex<double> > coldRoots( 
				double w, std::complex<double> ky, std::complex<double> kz, arma::cx_mat _epsilon)
{

		// see rsfxc_1D.nb for the derivation of these polynomial coeffs.

		std::complex<double> k4 = pow(w,2)/pow(_c,2) * _epsilon(0,0);

		std::complex<double> k3 = pow(w,2)/pow(_c,2) 
				* ( ky * ( _epsilon(0,1) + _epsilon(1,0) ) 
						+ kz * ( _epsilon(0,2) + _epsilon(2,0) ) );

		std::complex<double> k2 = pow(w,2)/pow(_c,4) * ( pow(_c,2) * pow(ky,2) * ( _epsilon(0,0) + _epsilon(1,1) ) 
						+ pow(_c,2) * ky * kz * ( _epsilon(1,2) + _epsilon(2,1) ) 
						+ pow(_c,2) * pow(kz,2) * ( _epsilon(0,0) + _epsilon(2,2) ) 
						+ pow(w,2) * ( _epsilon(0,1) * _epsilon(1,0) 
										+ _epsilon(0,2) * _epsilon(2,0)
										- _epsilon(0,0) * ( _epsilon(1,1) + _epsilon(2,2) ) ) );

		std::complex<double> k1 = pow(w,2)/pow(_c,4) * ( pow(_c,2) * pow(ky,3)*( _epsilon(0,1) + _epsilon(1,0) )
						+ pow(_c,2) * pow(ky,2)*kz* ( _epsilon(0,2) + _epsilon(2,0) )
						+ kz * ( pow(_c,2) * pow(kz,2) * ( _epsilon(0,2) + _epsilon(2,0) )
							+ pow(w,2) * ( _epsilon(0,1) * _epsilon(1,2) 
										- _epsilon(1,1) * ( _epsilon(0,2) + _epsilon(2,0) )
										+ _epsilon(1,0) * _epsilon(2,1) ) ) 
						+ ky * ( pow(w,2) * ( _epsilon(1,2) * _epsilon(2,0)
												+ _epsilon(0,2) * _epsilon(2,1) )
									+ ( _epsilon(0,1) + _epsilon(1,0) ) 
										* ( pow(_c,2) * pow(kz,2) - pow(w,2) * _epsilon(2,2) ) ) );
		
		std::complex<double> k0 = pow(w,2)/pow(_c,6) * ( pow(_c,4) * pow(ky,4) * _epsilon(1,1) 
							+ pow(_c,4) * kz * pow(ky,3) * ( _epsilon(1,2) + _epsilon(2,1) ) 
							+ pow(_c,2) * ky * kz 
                                * ( pow(_c,2) * pow(kz,2) * ( _epsilon(1,2) + _epsilon(2,1) ) 
									+ pow(w,2) * ( _epsilon(0,2) * _epsilon(1,0) 
													+ _epsilon(0,1) * _epsilon(2,0) 
													- _epsilon(0,0) * 
														( _epsilon(1,2) + _epsilon(2,1) ) ) ) 
							+ pow(_c,4) * pow(kz,4) * _epsilon(2,2) 
							+ pow(_c,2) * pow(w,2) * pow(kz,2) 
								* ( _epsilon(0,2) * _epsilon(2,0) 
									+ _epsilon(1,2) * _epsilon(2,1) 
									- ( _epsilon(0,0) + _epsilon(1,1) ) * _epsilon(2,2) ) 
							+pow(w,4) * ( _epsilon(0,2) * ( -_epsilon(1,1) * _epsilon(2,0) 
															+ _epsilon(1,0) * _epsilon(2,1) ) 
								+ _epsilon(0,1) * ( _epsilon(1,2) * _epsilon(2,0) 
													- _epsilon(1,0) * _epsilon(2,2) ) 
								+ _epsilon(0,0) * ( -_epsilon(1,2) * _epsilon(2,1) 
													+ _epsilon(1,1) * _epsilon(2,2) ) ) 
							+ pow(_c,2) * pow(ky,2) * ( pow(_c,2) * pow(kz,2) * 
								( _epsilon(1,1) + _epsilon(2,2) ) 
								+ pow(w,2) * ( _epsilon(0,1) * _epsilon(1,0) 
											+ _epsilon(1,2) * _epsilon(2,1) 
											- _epsilon(1,1) * 
												( _epsilon(0,0) + _epsilon(2,2) ) ) ) );
#ifdef _DEBUG
		std::cout << "Coeffs: " << std::endl;
		std::cout << "k0: " << k0 << std::endl;
		std::cout << "k1: " << k1 << std::endl;
		std::cout << "k2: " << k2 << std::endl;
		std::cout << "k3: " << k3 << std::endl;
		std::cout << "k4: " << k4 << std::endl;
#endif
		o2scl::simple_quartic_complex quart;
		std::vector<std::complex<double> > roots;
		roots.clear();
		roots.resize(4);
		std::complex<double> rootsTmp[4];
		quart.solve_c(k4,k3,k2,k1,k0,rootsTmp[0],rootsTmp[1],rootsTmp[2],rootsTmp[3]);

		for(int i=0;i<4;i++)
		{
				roots[i] = rootsTmp[i];
#ifdef _DEBUG
				std::cout << "Root: " << roots[i] << std::endl;
#endif
		}

		//std::sort(rootsVec.begin(),rootsVec.end());

		return(roots);
}

void DispersionApp::UpdateCalculation()
{
			//UpdateIonSpecies();
			UpdateBField();
			// Update number of plot data points

			std::cout << "Updating plots ..." << std::endl;

			int nRows = a_model->rowCount();
			if(nRows<nX) a_model->insertRows(nRows-1,nX-nRows);
			if(nRows>nX) a_model->removeRows(0,nRows-nX);
			std::cout << "nX: "<<nX<<" nRows: "<<a_model->rowCount()<<std::endl;

			if(nRows<nX) b_model->insertRows(nRows-1,nX-nRows);
			if(nRows>nX) b_model->removeRows(0,nRows-nX);

			for(int i=0;i<nX;i++)
			{

				float _bMag = bMag(i);
				std::vector<PlasmaSpecies> AllSpecies;
				std::vector<HotPlasmaSpecies> AllSpeciesHot;
				double _ionDensity = 0.0;
				//ions
				for(int s=0;s<nSpecies;s++)
				{
						double _z;
						try {
							_z = boost::lexical_cast<double>(IonSpecies[s].SpecZ_LineEdit->text().narrow());
							std::cout << "Success: SpecZ_LineEdit read as " << _z << std::endl;
						} catch(boost::bad_lexical_cast &) {
							std::cout << "Error: SpecZ_LineEdit could not be read" << std::endl;
						}

						double _amu;
						try {
					   		_amu = boost::lexical_cast<double>(IonSpecies[s].SpecAMU_LineEdit->text().narrow());
							std::cout << "Success: SpecAMU_LineEdit read as " << _amu << std::endl;
						} catch(boost::bad_lexical_cast &) {
							std::cout << "Error: SpecAMU_LineEdit could not be read" << std::endl;
						}

						double _T_eV;
						try {
					   		_T_eV = boost::lexical_cast<double>(IonSpecies[s].SpecT_eV_LineEdit->text().narrow());
							std::cout << "Success: SpecT_eV_LineEdit read as " << _T_eV << std::endl;
						} catch(boost::bad_lexical_cast &) {
							std::cout << "Error: SpecT_eV_LineEdit could not be read" << std::endl;
						}

						double _n;
						try {
					   		_n = boost::lexical_cast<double>(IonSpecies[s].SpecDensity_m3_LineEdit->text().narrow());
							std::cout << "Success: Density_m3_LineEdit read as " << _n << std::endl;
						} catch(boost::bad_lexical_cast &) {
							std::cout << "Error: Density_m3_LineEdit could not be read" << std::endl;
						}

						int _maxHarmN;
						try {
					   		_maxHarmN = boost::lexical_cast<int>(IonSpecies[s].SpecHarmonicN_LineEdit->text().narrow());
							std::cout << "Success: SpecHarmonicN_LineEdit read as " << _maxHarmN << std::endl;
						} catch(boost::bad_lexical_cast &) {
							std::cout << "Error: SpecHarmonicN_LineEdit could not be read" << std::endl;
						}


						AllSpecies.push_back(PlasmaSpecies(_z,_amu,_n,_bMag));
						AllSpeciesHot.push_back(HotPlasmaSpecies(_z,_amu,_n,_bMag,_T_eV,_maxHarmN));
						_ionDensity += _z * _n;
						std::cout << "_ionDensity: " << _ionDensity << std::endl;

						std::cout << "AllSpecies[" << s << "] wp = " << AllSpecies[s].wp << std::endl;
						std::cout << "AllSpecies[" << s << "] wc = " << AllSpecies[s].wc << std::endl;

				}
				//electrons
				double _T_eV_e;
				try {
					_T_eV_e = boost::lexical_cast<double>(electron_T_eV_LineEdit->text().narrow());
					std::cout << "Success: electron_T_eV_LineEdit read as " << _T_eV_e << std::endl;
				} catch(boost::bad_lexical_cast &) {
					std::cout << "Error: electron_T_eV_LineEdit could not be read" << std::endl;
				}

				double _amu_e;
				try {
					_amu_e = boost::lexical_cast<double>(electron_AMU_LineEdit->text().narrow());
					std::cout << "Success: electron_AMU_LineEdit read as " << _amu_e << std::endl;
				} catch(boost::bad_lexical_cast &) {
					std::cout << "Error: electron_AMU_LineEdit could not be read" << std::endl;
				}

				int _maxHarmN;
				try {
					_maxHarmN = boost::lexical_cast<int>(electron_HarmonicN_LineEdit->text().narrow());
					std::cout << "Success: electron_LineEdit read as " << _maxHarmN << std::endl;
				} catch(boost::bad_lexical_cast &) {
					std::cout << "Error: electron_LineEdit could not be read" << std::endl;
				}

				AllSpecies.push_back(PlasmaSpecies(-1.0,_amu_e,_ionDensity,_bMag));
				AllSpeciesHot.push_back(HotPlasmaSpecies(-1.0,_amu_e,_ionDensity,_bMag,_T_eV_e,_maxHarmN));

				std::cout << "AllSpecies[e] wp = " << AllSpecies[nSpecies].wp << std::endl;
				std::cout << "AllSpecies[e] wc = " << AllSpecies[nSpecies].wc << std::endl;

				//kPar = boost::lexical_cast<double>(kPar_LineEdit->text().narrow());
				kyRe = boost::lexical_cast<double>(kyRe_LineEdit->text().narrow());
				kzRe = boost::lexical_cast<double>(kzRe_LineEdit->text().narrow());
				kyIm = boost::lexical_cast<double>(kyIm_LineEdit->text().narrow());
				kzIm = boost::lexical_cast<double>(kzIm_LineEdit->text().narrow());

				ky = std::complex<double>(kyRe,kyIm);
				kz = std::complex<double>(kzRe,kzIm);

				double _freq = boost::lexical_cast<double>(Freq_Hz_LineEdit->text().narrow());
				double _nuOmg = boost::lexical_cast<double>(nuOmg_LineEdit->text().narrow());

				double _omega = 2.0 * _pi * _freq;
				std::complex<double> _omega_c(_omega,_omega*_nuOmg);
				StixVars stix(_omega_c, AllSpecies);

				std::cout << "Stix.R = " << stix.R << std::endl;
				std::cout << "Stix.L = " << stix.L << std::endl;
				std::cout << "Stix.D = " << stix.D << std::endl;
				std::cout << "Stix.S = " << stix.S << std::endl;
				std::cout << "Stix.P = " << stix.P << std::endl;

				arma::colvec b_xyz = arma::zeros<arma::colvec>(3);
				b_xyz(0) = bx(i);
				b_xyz(1) = by(i);
				b_xyz(2) = bz(i); // remember the z direction is -phi

				std::cout<<"Initializing rotation matricies ..."<<std::endl;
				RotationMatrix rot(b_xyz);

				std::cout<<"Initializing cold dielectric in stx frame ..."<<std::endl;
				dielectric epsilonCold(stix);
				epsilonCold.epsilon.print("Cold dielectric stx: ");

				// I think these rotations need to change, i.e., need
				// to make it go abp->stx->xyz, such that the cold rotaion
				// does NOT depend on k.

				//epsilonCold.rotate(rot.stx2abp);
				//epsilonCold.epsilon.print("Cold dielectric abp: ");

				epsilonCold.rotate(rot.abp2xyz);
				epsilonCold.epsilon.print("Cold dielectric xyz: ");

				std::vector<std::complex<double> > ColdRoots = coldRoots(_omega,ky,kz,epsilonCold.epsilon);
				for(int c=0;c<4;c++)
				{
						std::cout<<"cold roots: "<<ColdRoots[c]<<std::endl;
				}
	
				dielectric epsilonHot(AllSpeciesHot,_omega_c,3,ky,kz,rot);

				std::complex<double> kx_guess;
				kx_guess = std::complex<double>(9.77,0.0);

				epsilonHot.populateSwansonKs(kx_guess);
				epsilonHot.epsilon.print("Hot dielectric abp: ");
				epsilonHot.rotate(rot.abp2xyz);
				epsilonHot.epsilon.print("Hot dielectric xyz: ");

				std::vector<std::complex<double> > HotRoots = coldRoots(_omega,ky,kz,epsilonHot.epsilon);
				for(int c=0;c<4;c++)
				{
						std::cout<<"hot roots: "<<HotRoots[c]<<std::endl;
				}
	
				if(i==nX-1)
				{
		
					std::cout<<"i: "<<i<<std::endl;

					// Just map the freaking space out and contour it. Then just look at a single 
					// point at a time in the complex k space for zeros. 
					
					double kx_reMin;
					try {
						kx_reMin = boost::lexical_cast<double>(kxReMin_LineEdit->text().narrow());
						std::cout << "Success: kxReMin_LineEdit read as " << kx_reMin << std::endl;
					} catch(boost::bad_lexical_cast &) {
						std::cout << "Error: kxReMin_LineEdit could not be read" << std::endl;
					}

					double kx_reMax;
					try {
						kx_reMax = boost::lexical_cast<double>(kxReMax_LineEdit->text().narrow());
						std::cout << "Success: kxReMax_LineEdit read as " << kx_reMax << std::endl;
					} catch(boost::bad_lexical_cast &) {
						std::cout << "Error: kxReMax_LineEdit could not be read" << std::endl;
					}

					double kx_imMin;
					try {
						kx_imMin = boost::lexical_cast<double>(kxImMin_LineEdit->text().narrow());
						std::cout << "Success: kxImMin_LineEdit read as " << kx_imMin << std::endl;
					} catch(boost::bad_lexical_cast &) {
						std::cout << "Error: kxImMin_LineEdit could not be read" << std::endl;
					}

					double kx_imMax;
					try {
						kx_imMax = boost::lexical_cast<double>(kxImMax_LineEdit->text().narrow());
						std::cout << "Success: kxImMax_LineEdit read as " << kx_imMax << std::endl;
					} catch(boost::bad_lexical_cast &) {
						std::cout << "Error: kxImMax_LineEdit could not be read" << std::endl;
					}

					double _nX;
					try {
						_nX = boost::lexical_cast<int>(kxReN_LineEdit->text().narrow());
						std::cout << "Success: kxReN_LineEdit read as " << _nX << std::endl;
					} catch(boost::bad_lexical_cast &) {
						std::cout << "Error: kxReN_LineEdit could not be read" << std::endl;
					}

					double _nY;
					try {
						_nY = boost::lexical_cast<int>(kxImN_LineEdit->text().narrow());
						std::cout << "Success: kxImN_LineEdit read as " << _nY << std::endl;
					} catch(boost::bad_lexical_cast &) {
						std::cout << "Error: kxImN_LineEdit could not be read" << std::endl;
					}

					blitz::Array<blitz::TinyVector<double,3>, 2> arr(_nX,_nY);

					for(int ii=0;ii<_nX;ii++)
					{
						for(int jj=0;jj<_nY;jj++)
						{
							double kx_re = (kx_reMax-kx_reMin)/(_nX-1)*ii+kx_reMin;
							double kx_im = (kx_imMax-kx_imMin)/(_nY-1)*jj+kx_imMin;

							std::complex<double> kxIn = std::complex<double>(kx_re,kx_im);	
							//std::cout<<"kxIn: "<<kxIn<<std::endl;
							std::complex<double> det = epsilonHot.determinant(kxIn);

							arr(ii,jj)[0] = kx_re;
							arr(ii,jj)[1] = kx_im;
							arr(ii,jj)[2] = log10(std::abs(det));
							//std::cout<<arr(ii,jj)[2]<<std::endl;
							if(arr(ii,jj)[2]!=arr(ii,jj)[2])
							{
									std::cout<<"Error: Nans in determinant."<<std::endl;
									exit(1);
							}
						}
					}

					Gnuplot gp("gnuplot -persist");
					//Gnuplot gp(fopen("external_text.gnu","w"));
					gp<<"set contour base"<<std::endl;
					gp<<"unset clabel"<<std::endl;
					gp<<"set xlabel 'Re(kx)'"<<std::endl;
					gp<<"set ylabel 'Im(kx)'"<<std::endl;
					gp<<"set cntrparam linear"<<std::endl;
					gp<<"set cntrparam levels incremental 0.5,0.5,200"<<std::endl;
					gp<<"set key off"<<std::endl;
					gp<<"unset surface; set view map"<<std::endl;
					gp<<"splot '-' with lines lc '#000000'"<<std::endl;
					gp.send(arr);
					//gp<<"splot"<<gp.file(arr,"external_text.dat")<<"with lines lc '#000000'"<<std::endl;
					
					blitz::Array<blitz::TinyVector<double,3>, 2> arrCold(_nX,_nY);

					for(int ii=0;ii<_nX;ii++)
					{
						for(int jj=0;jj<_nY;jj++)
						{
							double kx_re = (kx_reMax-kx_reMin)/(_nX-1)*ii+kx_reMin;
							double kx_im = (kx_imMax-kx_imMin)/(_nY-1)*jj+kx_imMin;

							std::complex<double> kxIn = std::complex<double>(kx_re,kx_im);	
							int cold=1;
							std::complex<double> det = epsilonCold.determinant(kxIn,ky,kz,rot,_omega,cold);
							arrCold(ii,jj)[0] = kx_re;
							arrCold(ii,jj)[1] = kx_im;
							arrCold(ii,jj)[2] = log10(std::abs(det));
							if(arrCold(ii,jj)[2]!=arrCold(ii,jj)[2])
							{
									std::cout<<"Error: Nans in determinant."<<std::endl;
									exit(1);
							}
						}
					}

					Gnuplot gpCold("gnuplot -persist");
					//Gnuplot gpCold(fopen("external_text.gnu","w"));
					gpCold<<"set contour base"<<std::endl;
					gpCold<<"unset clabel"<<std::endl;
					gpCold<<"set xlabel 'Re(kx)'"<<std::endl;
					gpCold<<"set ylabel 'Im(kx)'"<<std::endl;
					gpCold<<"set cntrparam linear"<<std::endl;
					gpCold<<"set cntrparam levels incremental 0.5,0.5,200"<<std::endl;
					gpCold<<"set key off"<<std::endl;
					gpCold<<"unset surface; set view map"<<std::endl;
					gpCold<<"splot '-' with lines lc '#000000'"<<std::endl;
					gpCold.send(arrCold);
					//gpCold<<"splot"<<gpCold.file(arrCold,"external_text.dat")<<"with lines lc '#000000'"<<std::endl;

				}
				
				// add new points to plot
				//for (unsigned i = 0; i < nX; ++i) {
				a_model->setData(i, 0, x[i]);

				a_model->setData(i, 1, std::real(ColdRoots[0]));
		  		a_model->setData(i, 2, std::real(ColdRoots[1]));
				a_model->setData(i, 3, std::real(ColdRoots[2]));
		  		a_model->setData(i, 4, std::real(ColdRoots[3]));

				a_model->setData(i, 5, std::imag(ColdRoots[0]));
		  		a_model->setData(i, 6, std::imag(ColdRoots[1]));
				a_model->setData(i, 7, std::imag(ColdRoots[2]));
		  		a_model->setData(i, 8, std::imag(ColdRoots[3]));

				//}

				//for (unsigned i = 0; i < nX; ++i) {
				  b_model->setData(i, 0, x[i]);
				  b_model->setData(i, 1, bMag(i));
				//}

				//AllSpecies.clear();
				//AllSpeciesHot.clear();
			}

}



Wt::WApplication *createApplication(const Wt::WEnvironment& env)
{
		    return new DispersionApp(env);
}

int main(int argc, char **argv)
{
		    return Wt::WRun(argc, argv, &createApplication);
}
