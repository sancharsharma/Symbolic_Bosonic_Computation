
import sympy as sym
import sympy.physics.secondquant as SQ

class Des(SQ.AnnihilateBoson):

	op_symbol = 'Des'

	def _dagger_(self):
		return Cre(self.state)
	
	def __repr__(self):
		return "Des(%s)" % self.state

	def _latex(self,printer):
		if self.state == 0:
			return "\hat{b}_0"
		else:
			return "\hat{%s}" % self.state.name
	
	def __eq__(self,other):
		if isinstance(other,SQ.B):
			return self.state == other.state
		else:
			return False
	
	def __hash__(self):
		return hash('AnnihilationOperator' + self.state.name)


class Cre(SQ.CreateBoson):

	op_symbol = 'Cre'

	def _dagger_(self):
		return Des(self.state)

	def __repr__(self):
		return "Cre(%s)" % self.state

	def _latex(self,printer):
		if self.state == 0:
			return "{{\hat{b}^{\\dagger}}_0}"
		else:
			return "{\hat{%s}^{\\dagger}}" % self.state.name

	def __eq__(self,other):
		if isinstance(other,SQ.Bd):
			return self.state == other.state
		else:
			return False
	
	def __hash__(self):
		return hash('CreationOperator' + self.state.name)


class Num(SQ.SqOperator):
	
	op_symbol = 'Num'

	def _dagger_(self):
		return self

	def __repr__(self):
		return "Num(%s)" % self.state

	def _latex(self,printer):
		if self.state == 0:
			return "\hat{n}_0"
		else:
			return "\hat{n}_{%s}" % self.state.name

	def __eq__(self,other):
		if isinstance(other,self.__class__):
			return self.state.name == other.state.name
		else:
			return False
	
	def __hash__(self):
		return hash('NumberBosonic' + self.state.name)


class SQ_Symbol(sym.Symbol):	
	
	def _dagger_(self):
		return sym.conjugate(self)
