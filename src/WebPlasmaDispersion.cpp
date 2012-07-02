#include <Wt/WApplication>
#include <Wt/WBreak>
#include <Wt/WContainerWidget>
#include <Wt/WLineEdit>
#include <Wt/WPushButton>
#include <Wt/WText>
#include <Wt/WDoubleValidator>
#include <Wt/WValidator>
#include <Wt/WColor>

class HelloApplication : public Wt::WApplication
{
		public:
		    HelloApplication(const Wt::WEnvironment& env);

		private:
			Wt::WLineEdit *nameEdit_;
			Wt::WText *greeting_;
			Wt::WLineEdit *Density_m3_LineEdit, *BField_LineEdit;
			Wt::WDoubleValidator *DensityValidator, *BFieldValidator;
			Wt::WText *DensityValidText, *BFieldValidText;

			void greet();
			void editedLineEdit(Wt::WLineEdit *LineEditIn, Wt::WText *ValidText);
};

HelloApplication::HelloApplication(const Wt::WEnvironment& env) : Wt::WApplication(env)
{
		    setTitle("Plasma Dispersion Calculator");

			root()->addWidget(new Wt::WText("Your name, please ? "));
			nameEdit_ = new Wt::WLineEdit(root());
			Wt::WPushButton *button = new Wt::WPushButton("RUN", root());
			root()->addWidget(new Wt::WBreak());
			greeting_ = new Wt::WText(root());
			button->clicked().connect(this, &HelloApplication::greet);

			Wt::WText *DensityText = new Wt::WText(root());
			DensityText->setText("Density [1/m^3]");
			DensityValidator = new Wt::WDoubleValidator(0,5,root());
			DensityValidator->setMandatory(1);
			Density_m3_LineEdit = new Wt::WLineEdit(root());
			Density_m3_LineEdit->setValidator(DensityValidator);
			DensityValidText = new Wt::WText(root());
			Density_m3_LineEdit->changed().connect(
					boost::bind(&HelloApplication::editedLineEdit,this,Density_m3_LineEdit,DensityValidText));
			root()->addWidget(new Wt::WBreak());

			Wt::WText *BFieldText = new Wt::WText(root());
			BFieldText->setText("BField [T]");
			BFieldValidator = new Wt::WDoubleValidator(0,5,root());
			BFieldValidator->setMandatory(1);
			BField_LineEdit = new Wt::WLineEdit(root());
			BField_LineEdit->setValidator(BFieldValidator);
			BFieldValidText = new Wt::WText(root());
			BField_LineEdit->changed().connect(
					boost::bind(&HelloApplication::editedLineEdit,this,BField_LineEdit,BFieldValidText));

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
