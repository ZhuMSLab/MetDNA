import { MetDNAAngularPage } from './app.po';

describe('met-dna-angular App', () => {
  let page: MetDNAAngularPage;

  beforeEach(() => {
    page = new MetDNAAngularPage();
  });

  it('should display welcome message', () => {
    page.navigateTo();
    expect(page.getParagraphText()).toEqual('Welcome to app!');
  });
});
