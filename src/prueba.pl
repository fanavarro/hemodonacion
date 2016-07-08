#!/usr/bin/perl -w -- 
#
# generated by wxGlade 0.6.8 on Thu Jul  7 19:03:34 2016
#
# To get wxPerl visit http://wxPerl.sourceforge.net/
#

use Wx 0.15 qw[:allclasses];
use strict;

# begin wxGlade: dependencies
# end wxGlade

# begin wxGlade: extracode
# end wxGlade

package MainFrame;

use Wx qw[:everything];
use base qw(Wx::Frame);
use strict;

use Wx::Locale gettext => '_T';
sub new {
    my( $self, $parent, $id, $title, $pos, $size, $style, $name ) = @_;
    $parent = undef              unless defined $parent;
    $id     = -1                 unless defined $id;
    $title  = ""                 unless defined $title;
    $pos    = wxDefaultPosition  unless defined $pos;
    $size   = wxDefaultSize      unless defined $size;
    $name   = ""                 unless defined $name;

    # begin wxGlade: MainFrame::new
    $style = wxDEFAULT_FRAME_STYLE 
        unless defined $style;

    $self = $self->SUPER::new( $parent, $id, $title, $pos, $size, $style, $name );
    $self->{text_ctrl_1} = Wx::TextCtrl->new($self, wxID_ANY, "", wxDefaultPosition, wxDefaultSize, );
    $self->{run} = Wx::Button->new($self, wxID_ANY, _T("Go!"));
    $self->{result_textbox} = Wx::StaticText->new($self, wxID_ANY, _T("Results"), wxDefaultPosition, wxDefaultSize, );
    Wx::Event::EVT_BUTTON($b, -1, sub {
      my ($b, $evt) = @_;
      $b->SetLabel('Clicked');
      $b->Disable;
    });
    $self->__set_properties();
    $self->__do_layout();

    # end wxGlade
    return $self;

}


sub __set_properties {
    my $self = shift;
    # begin wxGlade: MainFrame::__set_properties
    $self->SetTitle(_T("Kozak Info"));
    $self->SetSize(Wx::Size->new(835, 749));
    # end wxGlade
}

sub __do_layout {
    my $self = shift;
    # begin wxGlade: MainFrame::__do_layout
    $self->{sizer_3} = Wx::BoxSizer->new(wxVERTICAL);
    $self->{sizer_5} = Wx::BoxSizer->new(wxHORIZONTAL);
    $self->{sizer_4} = Wx::BoxSizer->new(wxHORIZONTAL);
    $self->{sizer_4}->Add($self->{text_ctrl_1}, 0, 0, 0);
    $self->{sizer_4}->Add($self->{run}, 0, 0, 0);
    $self->{sizer_3}->Add($self->{sizer_4}, 1, wxEXPAND, 0);
    $self->{sizer_5}->Add($self->{result_textbox}, 0, wxEXPAND, 0);
    $self->{sizer_3}->Add($self->{sizer_5}, 1, wxEXPAND, 0);
    $self->SetSizer($self->{sizer_3});
    $self->Layout();
    # end wxGlade
}

# end of class MainFrame

1;

1;

package main;

unless(caller){
    my $local = Wx::Locale->new("Spanish", "es", "es"); # replace with ??
    $local->AddCatalog("app"); # replace with the appropriate catalog name

    local *Wx::App::OnInit = sub{1};
    my $app = Wx::App->new();
    Wx::InitAllImageHandlers();

    my $MainFrame = MainFrame->new();

    $app->SetTopWindow($MainFrame);
    $MainFrame->Show(1);
    $app->MainLoop();
}
