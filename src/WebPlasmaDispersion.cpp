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
#include "zcircle.hpp"

class IonSpec
{
		public:
		IonSpec(int, Wt::WGroupBox*);
		~IonSpec(void);
		void setVisible(void);
		void setHidden(void);

		private:

		public:
		Wt::WLineEdit *SpecAMU_LineEdit, *SpecZ_LineEdit;
		Wt::WBreak *SpecBreak;
		Wt::WText *SpecText;
		int SpecNo;
};

IonSpec::IonSpec(int SpecNoIn, Wt::WGroupBox *IonSpecies_Container )
{
	SpecNo = SpecNoIn;
	std::stringstream TmpStr;
	TmpStr << "Ion Species " << SpecNoIn << " (AMU, Z)";
	Wt::WString TmpStrWt(TmpStr.str());
	SpecText = new Wt::WText(TmpStrWt,IonSpecies_Container);
	SpecText->setHidden(1);
	SpecAMU_LineEdit = new Wt::WLineEdit("1",IonSpecies_Container);
	SpecAMU_LineEdit->setHidden(1);
	SpecZ_LineEdit = new Wt::WLineEdit("1",IonSpecies_Container);
	SpecZ_LineEdit->setHidden(1);
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
	SpecBreak->setHidden(0);
}
void IonSpec::setHidden(void)
{
	SpecText->setHidden(1);
	SpecAMU_LineEdit->setHidden(1);
	SpecZ_LineEdit->setHidden(1);
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
			Wt::WLineEdit *Density_m3_LineEdit, *BField_LineEdit, *Freq_Hz_LineEdit;
			Wt::WDoubleValidator *DensityValidator, *BFieldValidator, *FreqValidator;
			Wt::WText *DensityValidText, *BFieldValidText, *FreqValidText;

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
			float r0, b0, rMin, rMax, kPar, ky, kz, yRangeMin, yRangeMax;
			double Temp_eV;
			Wt::WText *r0_Text, *b0_Text, *rMin_Text, *rMax_Text, *kPar_Text,
					*ky_Text, *kz_Text, *Temp_eV_Text;
			Wt::WLineEdit *r0_LineEdit, *b0_LineEdit, *rMin_LineEdit, *rMax_LineEdit, 
					*kPar_LineEdit, *kz_LineEdit, *ky_LineEdit, *Temp_eV_LineEdit;

			// Plot Y-Range controllers
			Wt::WText *yRangeMin_Text, *yRangeMax_Text;
			Wt::WLineEdit *yRangeMin_LineEdit, *yRangeMax_LineEdit;

			std::vector <float> x,b;

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
			West_Container->resize(250,400);
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

			ParameterType_ButtonGroup->setCheckedButton(ParameterType_ButtonGroup->button(MagneticField));
			ParameterType_ButtonGroup->checkedChanged().connect(
							boost::bind(&DispersionApp::RevealParameterSettings,this));

			PointParameterType_Container = new Wt::WGroupBox("Point Parameter Settings",root());

			new Wt::WBreak(West_Container);
			Wt::WText *DensityText = new Wt::WText(West_Container);
			DensityText->setText("Density [1/m^3]");
			DensityValidator = new Wt::WDoubleValidator(1e1,1e25,West_Container);
			DensityValidator->setMandatory(1);
			Density_m3_LineEdit = new Wt::WLineEdit("2e19",West_Container);
			Density_m3_LineEdit->setValidator(DensityValidator);
			DensityValidText = new Wt::WText(West_Container);
			Density_m3_LineEdit->changed().connect(
					boost::bind(&DispersionApp::editedLineEdit,this,Density_m3_LineEdit,DensityValidText));

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

			PointParameterType_Container->setHidden(1);

			MagneticParameterType_Container = new Wt::WGroupBox("Magnetic Field Scan Settings",root());
			MagneticParameterType_Container->setHidden(0);

			//kPar_Text = new Wt::WText("kPar",West_Container);
			//kPar_LineEdit = new Wt::WLineEdit("0.1",West_Container);
			//new Wt::WBreak(West_Container);

			ky_Text = new Wt::WText("ky",West_Container);
			ky_LineEdit = new Wt::WLineEdit("0",West_Container);
			new Wt::WBreak(West_Container);

			kz_Text = new Wt::WText("kz",West_Container);
			kz_LineEdit = new Wt::WLineEdit("10",West_Container);
			new Wt::WBreak(West_Container);

			Temp_eV_Text = new Wt::WText("Temp_eV",West_Container);
			Temp_eV_LineEdit = new Wt::WLineEdit("0.01",West_Container);
			new Wt::WBreak(West_Container);

			r0_Text = new Wt::WText("r0",West_Container);
			r0_LineEdit = new Wt::WLineEdit("2",West_Container);
			new Wt::WBreak(West_Container);

			b0_Text = new Wt::WText("b0",West_Container);
			b0_LineEdit = new Wt::WLineEdit("5",West_Container);
			new Wt::WBreak(West_Container);

			rMin_Text = new Wt::WText("rMin",West_Container);
			rMin_LineEdit = new Wt::WLineEdit("0.1",West_Container);
			new Wt::WBreak(West_Container);

			rMax_Text = new Wt::WText("rMax",West_Container);
			rMax_LineEdit = new Wt::WLineEdit("4",West_Container);

			Update_PushButton = new Wt::WPushButton("Update",West_Container);
			Update_PushButton->clicked().connect(this, &DispersionApp::UpdateCalculation);

			IonSpecies_Container = new Wt::WGroupBox("Ion Species",root());

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

			nSpecies = 0;

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

			UpdateCalculation();
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
	_Chart->resize(900,600);
	_Chart->setPlotAreaPadding(200,Wt::Left);
	_Chart->setPlotAreaPadding(80,Wt::Bottom);
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
			r0 = boost::lexical_cast<float>(r0_LineEdit->text().narrow());

			x.resize(nX);
			b.resize(nX);

			for(int i=0;i<nX;i++)
			{
					x[i] = i*(rMax-rMin)/nX+rMin;
					b[i] = b0*r0/x[i];
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
				double w, std::complex<double> kz, std::complex<double> ky, arma::cx_mat _epsilon)
{

		// see rsfxc_1D.nb for the derivation of these polynomial coeffs.

		std::complex<double> k4 = pow(w,2)/pow(_c,2) * _epsilon(0,0);

		std::complex<double> k3 = pow(w,2)/pow(_c,2) 
				* ( (-kz) * ( _epsilon(0,1) + _epsilon(1,0) ) 
						+ ky * ( _epsilon(0,2) + _epsilon(2,0) ) );

		std::complex<double> k2 = pow(w,2)/pow(_c,4) * ( pow(_c,2) * pow((-kz),2) * ( _epsilon(0,0) + _epsilon(1,1) ) 
						+ pow(_c,2) * ky * (-kz) * ( _epsilon(1,2) + _epsilon(2,1) ) 
						+ pow(_c,2) * pow(ky,2) * ( _epsilon(0,0) + _epsilon(2,2) ) 
						+ pow(w,2) * ( _epsilon(0,1) * _epsilon(1,0) 
										+ _epsilon(0,2) * _epsilon(2,0)
										- _epsilon(0,0) * ( _epsilon(1,1) + _epsilon(2,2) ) ) );

		std::complex<double> k1 = pow(w,2)/pow(_c,4) * ( pow(_c,2) * pow(ky,2) * (-kz) * ( _epsilon(0,1) + _epsilon(1,0) )
						+ pow(_c,2) * pow(ky,3) * ( _epsilon(0,2) + _epsilon(2,0) )
						+ ky * ( pow(_c,2) * pow((-kz),2) * ( _epsilon(0,2) + _epsilon(2,0) )
							+ pow(w,2) * ( _epsilon(0,1) * _epsilon(1,2) 
										- _epsilon(1,1) * ( _epsilon(0,2) + _epsilon(2,0) )
										+ _epsilon(1,0) * _epsilon(2,1) ) ) 
						+ (-kz) * ( pow(w,2) * ( _epsilon(1,2) * _epsilon(2,0)
												+ _epsilon(0,2) * _epsilon(2,1) )
									+ ( _epsilon(0,1) + _epsilon(1,0) ) 
										* ( pow(_c,2) * pow((-kz),2) - pow(w,2) * _epsilon(2,2) ) ) );
		
		std::complex<double> k0 = pow(w,2)/pow(_c,6) * ( pow(_c,4) * pow((-kz),4) * _epsilon(1,1) 
							+ pow(_c,4) * ky * pow((-kz),3) * ( _epsilon(1,2) + _epsilon(2,1) ) 
							+ pow(_c,2) * ky * (-kz) 
                                * ( pow(_c,2) * pow(ky,2) * ( _epsilon(1,2) + _epsilon(2,1) ) 
									+ pow(w,2) * ( _epsilon(0,2) * _epsilon(1,0) 
													+ _epsilon(0,1) * _epsilon(2,0) 
													- _epsilon(0,0) * 
														( _epsilon(1,2) + _epsilon(2,1) ) ) ) 
							+ pow(_c,4) * pow(ky,4) * _epsilon(2,2) 
							+ pow(_c,2) * pow(w,2) * pow(ky,2) 
								* ( _epsilon(0,2) * _epsilon(2,0) 
									+ _epsilon(1,2) * _epsilon(2,1) 
									- ( _epsilon(0,0) + _epsilon(1,1) ) * _epsilon(2,2) ) 
							+pow(w,4) * ( _epsilon(0,2) * ( -_epsilon(1,1) * _epsilon(2,0) 
															+ _epsilon(1,0) * _epsilon(2,1) ) 
								+ _epsilon(0,1) * ( _epsilon(1,2) * _epsilon(2,0) 
													- _epsilon(1,0) * _epsilon(2,2) ) 
								+ _epsilon(0,0) * ( -_epsilon(1,2) * _epsilon(2,1) 
													+ _epsilon(1,1) * _epsilon(2,2) ) ) 
							+ pow(_c,2) * pow((-kz),2) * ( pow(_c,2) * pow(ky,2) * 
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
			UpdateIonSpecies();
			UpdateBField();
			// Update number of plot data points

			std::cout << "Updating plots ..." << std::endl;

			int nRows = a_model->rowCount();
			if(nRows<nX) a_model->insertRows(nRows-1,nX-nRows);
			if(nRows>nX) a_model->removeRows(0,nRows-nX);
			std::cout << "nX: "<<nX<<" nRows: "<<a_model->rowCount()<<std::endl;

			if(nRows<nX) b_model->insertRows(nRows-1,nX-nRows);
			if(nRows>nX) b_model->removeRows(0,nRows-nX);


			std::vector<double> tmp;

			for(int i=0;i<nX;i++)
			{

				float _bMag = b[i];
				std::vector<PlasmaSpecies> AllSpecies;
				std::vector<HotPlasmaSpecies> AllSpeciesHot;
				double _ionDensity = 0.0;
				Temp_eV = boost::lexical_cast<double>(Temp_eV_LineEdit->text().narrow());
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

						double _n;
						try {
					   		_n = boost::lexical_cast<double>(Density_m3_LineEdit->text().narrow());
							std::cout << "Success: Density_m3_LineEdit read as " << _n << std::endl;
						} catch(boost::bad_lexical_cast &) {
							std::cout << "Error: Density_m3_LineEdit could not be read" << std::endl;
						}

						AllSpecies.push_back(PlasmaSpecies(_z,_amu,_n,_bMag));
						AllSpeciesHot.push_back(HotPlasmaSpecies(_z,_amu,_n,_bMag,Temp_eV));
						_ionDensity += _z * _n;
						std::cout << "_ionDensity: " << _ionDensity << std::endl;

						std::cout << "AllSpecies[" << s << "] wp = " << AllSpecies[s].wp << std::endl;
						std::cout << "AllSpecies[" << s << "] wc = " << AllSpecies[s].wc << std::endl;

				}
				//electrons
				AllSpecies.push_back(PlasmaSpecies(-1.0,_me_mi,_ionDensity,_bMag));
				AllSpeciesHot.push_back(HotPlasmaSpecies(-1.0,_me_mi,_ionDensity,_bMag,Temp_eV));

				std::cout << "AllSpecies[e] wp = " << AllSpecies[nSpecies].wp << std::endl;
				std::cout << "AllSpecies[e] wc = " << AllSpecies[nSpecies].wc << std::endl;

				//kPar = boost::lexical_cast<double>(kPar_LineEdit->text().narrow());
				ky = boost::lexical_cast<double>(ky_LineEdit->text().narrow());
				kz = boost::lexical_cast<double>(kz_LineEdit->text().narrow());

				double _freq = boost::lexical_cast<double>(Freq_Hz_LineEdit->text().narrow());
				double _omega = 2.0 * _pi * _freq;
				StixVars stix(_omega, AllSpecies);

				std::cout << "Stix.R = " << stix.R << std::endl;
				std::cout << "Stix.L = " << stix.L << std::endl;
				std::cout << "Stix.D = " << stix.D << std::endl;
				std::cout << "Stix.S = " << stix.S << std::endl;
				std::cout << "Stix.P = " << stix.P << std::endl;

				tmp.push_back(stix.R);

				arma::colvec b_xyz = arma::zeros<arma::colvec>(3);
				b_xyz(0) = 0;
				b_xyz(1) = 0;
				b_xyz(2) = _bMag; // remember the z direction is -phi

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

				std::complex<double> ky_(ky,0.0);
				std::complex<double> kz_(kz,0.0);

				std::vector<std::complex<double> > ColdRoots = coldRoots(_omega,ky_,kz_,epsilonCold.epsilon);
				//std::vector<std::complex<double> > HotRoots = ColdRoots;
				
	
				dielectric epsilonHot(AllSpeciesHot,_omega,3,ky,kz,rot);

				std::complex<double> kx_guess;
				kx_guess = std::complex<double>(6.0,0.0);

				epsilonHot.populateSwansonKs(kx_guess);
				epsilonHot.rotate(rot.abp2xyz);
				epsilonHot.epsilon.print("Hot dielectric xyz: ");

				double w=25, e=0.9;
				int n=5, m=6;
				double root_re, root_im;

				std::cout<<"Calling Zcircle ..."<<std::endl;
				Zcircle(w,e,std::real(kx_guess),std::imag(kx_guess),n,m,root_re,root_im,epsilonHot);				

				//// starting from each cold root
				//for(int t=0;t<4;t++)
				//{
				//	arma::cx_colvec k_xyz(3);
				//	std::complex<double> krGuess = ColdRoots[t];

				//	std::cout<<std::endl<<std::endl;
				//	std::cout<<"original cold krGuess: "<<krGuess<<std::endl;
				//	// iterate to the hot root
				//	for(int it=0;it<40;it++)
				//	{
				//		k_xyz(0) = krGuess;
				//		k_xyz(1) = ky_;
				//		k_xyz(2) = kz_;

				//		//arma::cx_colvec k_abp = rot.xyz2abp * k_xyz;
				//		//rot.xyz2abp.print("xyz2abp: ");
				//		//k_xyz.print("k_xyz: ");
				//		//k_abp.print("k_abp: ");

				//		dielectric epsilonHot(AllSpeciesHot,_omega,3,ky,kz,rot);
				//		//epsilonHot.epsilon.print("Hot dielectric abp: ");

				//		epsilonHot.rotate(rot.abp2xyz);
				//		//epsilonHot.epsilon.print("Hot dielectric xyz: ");
				//		//epsilonCold.epsilon.print("Cold dielectric xyz: ");

				//		std::vector<std::complex<double> > HotRootsTmp = coldRoots(_omega,ky_,kz_,epsilonHot.epsilon);
				//		//std::cout<<"krIter: "<<HotRootsTmp[0]<<"   "<<ColdRoots[0]<<std::endl;
				//		//std::cout<<"krIter: "<<HotRootsTmp[1]<<"   "<<ColdRoots[1]<<std::endl;
				//		//std::cout<<"krIter: "<<HotRootsTmp[2]<<"   "<<ColdRoots[2]<<std::endl;
				//		//std::cout<<"krIter: "<<HotRootsTmp[3]<<"   "<<ColdRoots[3]<<std::endl;

				//		// find the same solution
				//		std::vector<double> diffVec(4,0);
				//		diffVec[0] = std::abs(std::imag(HotRootsTmp[0])-std::imag(krGuess));
				//		diffVec[1] = std::abs(std::imag(HotRootsTmp[1])-std::imag(krGuess));
				//		diffVec[2] = std::abs(std::imag(HotRootsTmp[2])-std::imag(krGuess));
				//		diffVec[3] = std::abs(std::imag(HotRootsTmp[3])-std::imag(krGuess));

				//		int idx = std::min_element(diffVec.begin(),diffVec.end())-diffVec.begin();

				//		arma::cx_colvec hotRoots = arma::zeros<arma::cx_colvec>(4);
				//		hotRoots(0) = HotRootsTmp[0];
				//		hotRoots(1) = HotRootsTmp[1];
				//		hotRoots(2) = HotRootsTmp[2];
				//		hotRoots(3) = HotRootsTmp[3];
				//		arma::colvec diff = arma::abs(arma::imag(hotRoots)-std::imag(krGuess));
				//		//hotRoots.print("hotRoots: ");
				//		//diff.print("diff: ");
				//		//std::cout<<"Min diff idx: "<<idx<<std::endl;
				//		krGuess = HotRootsTmp[idx];
				//		std::cout<<"New krGuess: "<<krGuess<<std::endl;

				//	}
				//}


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
				  b_model->setData(i, 1, b[i]);
				//}

				AllSpecies.clear();
				AllSpeciesHot.clear();
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
