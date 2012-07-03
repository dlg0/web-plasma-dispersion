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

class HelloApplication : public Wt::WApplication
{
		public:
		    HelloApplication(const Wt::WEnvironment& env);

		private:
			Wt::WLineEdit *nameEdit_;
			Wt::WText *greeting_;
			Wt::WLineEdit *Density_m3_LineEdit, *BField_LineEdit, *Freq_Hz_LineEdit;
			Wt::WDoubleValidator *DensityValidator, *BFieldValidator, *FreqValidator;
			Wt::WText *DensityValidText, *BFieldValidText, *FreqValidText;
			Wt::WLineEdit *Spec1AMU, *Spec2AMU, *Spec3AMU, *Spec4AMU, *Spec5AMU, *Spec6AMU;
			Wt::WLineEdit *Spec1Z, *Spec2Z, *Spec3Z, *Spec4Z, *Spec5Z, *Spec6Z;
			Wt::WText *Spec1Text, *Spec2Text, *Spec3Text, *Spec4Text, *Spec5Text, *Spec6Text;

			void greet();
			void editedLineEdit(Wt::WLineEdit *LineEditIn, Wt::WText *ValidText);
			void showIonSpecies(Wt::WComboBox *ComboBoxIn);

};

HelloApplication::HelloApplication(const Wt::WEnvironment& env) : Wt::WApplication(env)
{
		    setTitle("Plasma Dispersion Calculator");

			//root()->addWidget(new Wt::WText("Your name, please ? "));
			//nameEdit_ = new Wt::WLineEdit(root());
			//Wt::WPushButton *button = new Wt::WPushButton("RUN", root());
			//root()->addWidget(new Wt::WBreak());
			//greeting_ = new Wt::WText(root());
			//button->clicked().connect(this, &HelloApplication::greet);

			Wt::WText *DensityText = new Wt::WText(root());
			DensityText->setText("Density [1/m^3]");
			DensityValidator = new Wt::WDoubleValidator(1e1,1e25,root());
			DensityValidator->setMandatory(1);
			Density_m3_LineEdit = new Wt::WLineEdit("2e19",root());
			Density_m3_LineEdit->setValidator(DensityValidator);
			DensityValidText = new Wt::WText(root());
			Density_m3_LineEdit->changed().connect(
					boost::bind(&HelloApplication::editedLineEdit,this,Density_m3_LineEdit,DensityValidText));

			root()->addWidget(new Wt::WBreak());

			Wt::WText *BFieldText = new Wt::WText(root());
			BFieldText->setText("BField [T]");
			BFieldValidator = new Wt::WDoubleValidator(1e-5,5e5,root());
			BFieldValidator->setMandatory(1);
			BField_LineEdit = new Wt::WLineEdit("5",root());
			BField_LineEdit->setValidator(BFieldValidator);
			BFieldValidText = new Wt::WText(root());
			BField_LineEdit->changed().connect(
					boost::bind(&HelloApplication::editedLineEdit,this,BField_LineEdit,BFieldValidText));

			root()->addWidget(new Wt::WBreak());

			Wt::WText *FreqText = new Wt::WText(root());
			FreqText->setText("Freq [Hz]");
			FreqValidator = new Wt::WDoubleValidator(1e1,5e12,root());
			FreqValidator->setMandatory(1);
			Freq_Hz_LineEdit = new Wt::WLineEdit("30e6",root());
			Freq_Hz_LineEdit->setValidator(FreqValidator);
			FreqValidText = new Wt::WText(root());
			Freq_Hz_LineEdit->changed().connect(
					boost::bind(&HelloApplication::editedLineEdit,this,Freq_Hz_LineEdit,FreqValidText));

			root()->addWidget(new Wt::WBreak());

			Wt::WText *nIonSpecText = new Wt::WText("N Ion Species",root());
			Wt::WComboBox *nIonSpecies_ComboBox = new Wt::WComboBox(root());
			nIonSpecies_ComboBox->addItem("1");
			nIonSpecies_ComboBox->addItem("2");
			nIonSpecies_ComboBox->addItem("3");
			nIonSpecies_ComboBox->addItem("4");
			nIonSpecies_ComboBox->addItem("5");
			nIonSpecies_ComboBox->addItem("6");
			nIonSpecies_ComboBox->setCurrentIndex(0);
			nIonSpecies_ComboBox->changed().connect(
				boost::bind(&HelloApplication::showIonSpecies,this, nIonSpecies_ComboBox));

			root()->addWidget(new Wt::WBreak());

			Spec1Text = new Wt::WText("Ion Species 1 (AMU, Z)", root() );
			Spec1AMU = new Wt::WLineEdit("2",root());
			Spec1Z = new Wt::WLineEdit("1",root());
			root()->addWidget(new Wt::WBreak());

			Spec2Text = new Wt::WText("Ion Species 2 (AMU, Z)", root() );
			Spec2Text->setHidden(1);
			Spec2AMU = new Wt::WLineEdit(root());
			Spec2AMU->setHidden(1);
			Spec2Z = new Wt::WLineEdit(root());
			Spec2Z->setHidden(1);
			root()->addWidget(new Wt::WBreak());

			Spec3Text = new Wt::WText("Ion Species 3 (AMU, Z)", root() );
			Spec3Text->setHidden(1);
			Spec3AMU = new Wt::WLineEdit(root());
			Spec3AMU->setHidden(1);
			Spec3Z = new Wt::WLineEdit(root());
			Spec3Z->setHidden(1);
			root()->addWidget(new Wt::WBreak());

			Spec4Text = new Wt::WText("Ion Species 4 (AMU, Z)", root() );
			Spec4Text->setHidden(1);
			Spec4AMU = new Wt::WLineEdit(root());
			Spec4AMU->setHidden(1);
			Spec4Z = new Wt::WLineEdit(root());
			Spec4Z->setHidden(1);
			root()->addWidget(new Wt::WBreak());

			Spec5Text = new Wt::WText("Ion Species 5 (AMU, Z)", root() );
			Spec5Text->setHidden(1);
			Spec5AMU = new Wt::WLineEdit(root());
			Spec5AMU->setHidden(1);
			Spec5Z = new Wt::WLineEdit(root());
			Spec5Z->setHidden(1);
			root()->addWidget(new Wt::WBreak());

			Spec6Text = new Wt::WText("Ion Species 6 (AMU, Z)", root() );
			Spec6Text->setHidden(1);
			Spec6AMU = new Wt::WLineEdit(root());
			Spec6AMU->setHidden(1);
			Spec6Z = new Wt::WLineEdit(root());
			Spec6Z->setHidden(1);
			root()->addWidget(new Wt::WBreak());
}

void HelloApplication::showIonSpecies(Wt::WComboBox *ComboBoxIn)
{

		int nSpec = ComboBoxIn->currentIndex();
		std::cout << nSpec << std::endl;
		nSpec++;

		if(nSpec>=1){
			Spec1Text->setHidden(0);
			Spec1AMU->setHidden(0);
			Spec1Z->setHidden(0);
		}
		else
		{
			Spec1Text->setHidden(1);
			Spec1AMU->setHidden(1);
			Spec1Z->setHidden(1);
		}

		if(nSpec>=2){
			Spec2Text->setHidden(0);
			Spec2AMU->setHidden(0);
			Spec2Z->setHidden(0);
		}
		else
		{
			Spec2Text->setHidden(1);
			Spec2AMU->setHidden(1);
			Spec2Z->setHidden(1);
		}
	
		if(nSpec>=3){
			Spec3Text->setHidden(0);
			Spec3AMU->setHidden(0);
			Spec3Z->setHidden(0);
		}
		else
		{
			Spec3Text->setHidden(1);
			Spec3AMU->setHidden(1);
			Spec3Z->setHidden(1);
		}
	
		if(nSpec>=4){
			Spec4Text->setHidden(0);
			Spec4AMU->setHidden(0);
			Spec4Z->setHidden(0);
		}
		else
		{
			Spec4Text->setHidden(1);
			Spec4AMU->setHidden(1);
			Spec4Z->setHidden(1);
		}
		
		if(nSpec>=5){
			Spec5Text->setHidden(0);
			Spec5AMU->setHidden(0);
			Spec5Z->setHidden(0);
		}
		else
		{
			Spec5Text->setHidden(1);
			Spec5AMU->setHidden(1);
			Spec5Z->setHidden(1);
		}

		if(nSpec>=6){
			Spec6Text->setHidden(0);
			Spec6AMU->setHidden(0);
			Spec6Z->setHidden(0);
		}
		else
		{
			Spec6Text->setHidden(1);
			Spec6AMU->setHidden(1);
			Spec6Z->setHidden(1);
		}
	
}

void HelloApplication::editedLineEdit(Wt::WLineEdit *LineEditIn, Wt::WText *ValidText)
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

void HelloApplication::greet()
{
		    greeting_->setText("Hello there, " + nameEdit_->text());
}

Wt::WApplication *createApplication(const Wt::WEnvironment& env)
{
		    return new HelloApplication(env);
}

int main(int argc, char **argv)
{
		    return Wt::WRun(argc, argv, &createApplication);
}
